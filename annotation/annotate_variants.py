#!/usr/bin/env python3
"""
annotate_variants.py — Multi-database variant annotation orchestrator

Parses a VCF file (somatic or germline) and annotates each PASS variant
by chaining five annotation sources:

  1. VEP REST API   → gene, consequence, HGVS, SIFT/PolyPhen, regulatory
  2. gnomAD GraphQL → 8-population allele frequencies, popmax, faf95
  3. ClinVar eUtils → clinical significance, disease, review stars
  4. CADD web API   → deleteriousness score (raw + phred)
  5. CIViC GraphQL  → cancer drug evidence (somatic only)

Additionally normalizes HGVS notation via hgvs_formatter.py.

Usage:
    python annotation/annotate_variants.py \\
        --vcf results/somatic/filtered.vcf.gz \\
        --type somatic \\
        --output results/somatic/annotated.tsv \\
        --vcf-output results/somatic/annotated.vcf.gz

Output TSV columns (in order):
  COORDINATES:  chrom, pos, ref, alt
  GENE:         gene_symbol, ensembl_gene_id, gene_biotype, transcript_id, is_canonical, mane_select
  CONSEQUENCE:  consequence, impact, exon, intron
  HGVS:         hgvsg, hgvsc, hgvsp
  AMINO ACID:   amino_acids, codons, cdna_position, cds_position, protein_position
  IN SILICO:    sift_score, sift_pred, polyphen_score, polyphen_pred
  CONSERVATION: phylop, gerp, cadd_raw, cadd_phred
  POPULATION:   gnomad_af, gnomad_popmax_af, gnomad_popmax_pop,
                gnomad_afr_af, gnomad_amr_af, gnomad_asj_af, gnomad_eas_af,
                gnomad_fin_af, gnomad_nfe_af, gnomad_sas_af, gnomad_mid_af
  PATHOGENICITY:clinvar_clnsig, clinvar_review_stars, clinvar_clndn
  ACTIONABILITY:civic_evidence_level, civic_evidence_type, civic_drug, civic_cancer_type
  VCF FIELDS:   filter, vaf, depth, genotype (type-dependent)
  KNOWN:        existing_variation, regulatory_id
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

import cyvcf2

from vep_annotator    import VEPAnnotator
from gnomad_annotator import GnomADAnnotator
from clinvar_annotator import ClinVarAnnotator
from cadd_annotator   import CADDAnnotator
from civic_annotator  import CIViCAnnotator
from hgvs_formatter   import format_all, describe_variant

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Column order for output TSV ───────────────────────────────────────────────
COORDINATE_COLS = ["chrom", "pos", "ref", "alt", "filter"]
GENE_COLS       = ["gene_symbol", "ensembl_gene_id", "gene_biotype",
                   "transcript_id", "is_canonical", "mane_select"]
CONSEQUENCE_COLS= ["consequence", "impact", "exon", "intron"]
HGVS_COLS       = ["hgvsg", "hgvsc", "hgvsp"]
AA_COLS         = ["amino_acids", "codons", "cdna_position", "cds_position", "protein_position"]
INSILICO_COLS   = ["sift_score", "sift_pred", "polyphen_score", "polyphen_pred"]
CONSERV_COLS    = ["phylop", "gerp", "cadd_raw", "cadd_phred"]
POPULATION_COLS = ["gnomad_af", "gnomad_popmax_af", "gnomad_popmax_pop",
                   "gnomad_afr_af", "gnomad_amr_af", "gnomad_asj_af",
                   "gnomad_eas_af", "gnomad_fin_af", "gnomad_nfe_af",
                   "gnomad_sas_af", "gnomad_mid_af"]
CLINVAR_COLS    = ["clinvar_clnsig", "clinvar_review_stars", "clinvar_clndn",
                   "clinvar_review_status", "clinvar_last_evaluated"]
CIVIC_COLS      = ["civic_evidence_level", "civic_evidence_type",
                   "civic_clinical_significance", "civic_drug", "civic_cancer_type"]
VCF_COLS        = ["vaf", "depth", "genotype"]
KNOWN_COLS      = ["existing_variation", "regulatory_id", "af_1kg"]

ALL_COLS = (COORDINATE_COLS + GENE_COLS + CONSEQUENCE_COLS + HGVS_COLS
            + AA_COLS + INSILICO_COLS + CONSERV_COLS + POPULATION_COLS
            + CLINVAR_COLS + CIVIC_COLS + VCF_COLS + KNOWN_COLS)


def parse_vcf(vcf_path: str, variant_type: str) -> list[dict]:
    """
    Parse a VCF file into a list of variant dicts.
    Extracts FORMAT fields relevant to the variant type.
    """
    variants: list[dict] = []
    vcf = cyvcf2.VCF(vcf_path)

    for record in vcf:
        # Only process PASS variants
        if record.FILTER and record.FILTER != "PASS":
            continue

        base = {
            "chrom": record.CHROM,
            "pos":   record.POS,
            "ref":   record.REF,
            "alt":   record.ALT[0] if record.ALT else ".",
            "filter": record.FILTER or "PASS",
        }

        # Extract call-level fields
        if variant_type == "somatic":
            # For somatic: extract VAF from tumor sample (first sample with AF tag)
            try:
                af = record.format("AF")
                base["vaf"] = f"{float(af[0][0]):.4f}" if af is not None else "."
            except Exception:
                base["vaf"] = "."
            try:
                dp = record.format("DP")
                base["depth"] = str(int(dp[0][0])) if dp is not None else "."
            except Exception:
                base["depth"] = "."
            base["genotype"] = "somatic"

        else:  # germline
            # For germline: report genotypes for all trio members
            try:
                gts = [s["GT"] for s in record.samples]
                base["genotype"] = ",".join(
                    "/".join(str(a) for a in gt) for gt in gts
                )
            except Exception:
                base["genotype"] = "."
            try:
                dp = record.format("DP")
                depths = [str(int(d[0])) for d in dp] if dp is not None else ["."]
                base["depth"] = ",".join(depths)
            except Exception:
                base["depth"] = "."
            base["vaf"] = "."

        variants.append(base)

    vcf.close()
    logger.info("Parsed %d PASS variants from %s", len(variants), vcf_path)
    return variants


def write_tsv(variants: list[dict], output_path: str) -> None:
    """Write annotated variants to a tab-separated file."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_COLS, delimiter="\t",
                                extrasaction="ignore", restval=".")
        writer.writeheader()
        writer.writerows(variants)

    logger.info("Wrote %d variants to %s", len(variants), output_path)


def write_annotated_vcf(variants: list[dict], original_vcf: str, output_path: str) -> None:
    """
    Write an annotated VCF by adding INFO fields to the original VCF.
    INFO fields added: CSQ_GENE, CSQ_CONSEQUENCE, HGVSC, HGVSP, GNOMAD_AF,
                       CADD_PHRED, CLINVAR_CLNSIG, CIVIC_DRUG
    """
    info_fields = [
        ('CSQ_GENE',        'String',  'Gene symbol from VEP'),
        ('CSQ_CONSEQUENCE', 'String',  'VEP consequence (SO term)'),
        ('CSQ_IMPACT',      'String',  'VEP consequence impact (HIGH/MODERATE/LOW/MODIFIER)'),
        ('HGVSC',           'String',  'HGVS coding DNA notation'),
        ('HGVSP',           'String',  'HGVS protein notation'),
        ('HGVSG',           'String',  'HGVS genomic notation'),
        ('GNOMAD_AF',       'Float',   'gnomAD global allele frequency'),
        ('GNOMAD_POPMAX',   'Float',   'gnomAD popmax allele frequency'),
        ('CADD_PHRED',      'Float',   'CADD phred score'),
        ('CLINVAR_CLNSIG',  'String',  'ClinVar clinical significance'),
        ('CLINVAR_STARS',   'Integer', 'ClinVar review star rating (0-4)'),
        ('CIVIC_DRUG',      'String',  'CIViC associated drug (somatic)'),
    ]

    # vcf_hdr: used only to build the Writer header (Writer init may advance the handle)
    vcf_hdr = cyvcf2.VCF(original_vcf)
    for name, vtype, desc in info_fields:
        vcf_hdr.add_info_to_header({'ID': name, 'Number': 1, 'Type': vtype, 'Description': desc})

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    vcf_out = cyvcf2.Writer(output_path, vcf_hdr)
    vcf_hdr.close()

    # vcf_iter: fresh handle for iteration — must also have the new headers so
    # record.INFO[field] = value works (records inherit their VCF's header)
    vcf_iter = cyvcf2.VCF(original_vcf)
    for name, vtype, desc in info_fields:
        vcf_iter.add_info_to_header({'ID': name, 'Number': 1, 'Type': vtype, 'Description': desc})

    variant_map = {(v["chrom"], int(v["pos"]), v["ref"], v["alt"]): v for v in variants}

    for record in vcf_iter:
        key = (record.CHROM, record.POS, record.REF, record.ALT[0] if record.ALT else ".")
        ann = variant_map.get(key, {})

        def _f(v: str) -> float | None:
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        def _i(v: str) -> int | None:
            try:
                return int(v)
            except (TypeError, ValueError):
                return None

        str_fields = [
            ("CSQ_GENE",       ann.get("gene_symbol", ".")),
            ("CSQ_CONSEQUENCE",ann.get("consequence", ".")),
            ("CSQ_IMPACT",     ann.get("impact", ".")),
            ("HGVSC",          ann.get("hgvsc", ".")),
            ("HGVSP",          ann.get("hgvsp", ".")),
            ("HGVSG",          ann.get("hgvsg", ".")),
            ("CLINVAR_CLNSIG", ann.get("clinvar_clnsig", ".")),
            ("CIVIC_DRUG",     ann.get("civic_drug", ".")),
        ]
        for field, value in str_fields:
            if value and value != ".":
                record.INFO[field] = str(value)

        for field, raw, coerce in [
            ("GNOMAD_AF",    ann.get("gnomad_af", "."),         _f),
            ("GNOMAD_POPMAX",ann.get("gnomad_popmax_af", "."),  _f),
            ("CADD_PHRED",   ann.get("cadd_phred", "."),        _f),
            ("CLINVAR_STARS",ann.get("clinvar_review_stars","."),_i),
        ]:
            typed = coerce(raw) if raw and raw != "." else None
            if typed is not None:
                record.INFO[field] = typed

        vcf_out.write_record(record)

    vcf_out.close()
    vcf_iter.close()
    logger.info("Annotated VCF written to %s", output_path)


def save_checkpoint(variants: list[dict], path: str) -> None:
    with open(path, "w") as fh:
        json.dump(variants, fh)
    logger.info("  Checkpoint saved: %s", path)


def load_checkpoint(path: str) -> list[dict] | None:
    p = Path(path)
    if p.exists() and p.stat().st_size > 0:
        with open(path) as fh:
            data = json.load(fh)
        logger.info("  Resuming from checkpoint: %s (%d variants)", path, len(data))
        return data
    return None


def main():
    parser = argparse.ArgumentParser(description="Multi-database variant annotation pipeline")
    parser.add_argument("--vcf",        required=True,  help="Input VCF (PASS variants)")
    parser.add_argument("--type",       required=True,  choices=["somatic", "germline"])
    parser.add_argument("--output",     required=True,  help="Output TSV path")
    parser.add_argument("--vcf-output", default=None,   help="Output annotated VCF path")
    parser.add_argument("--threads",    type=int, default=1)
    parser.add_argument("--skip-cadd",  action="store_true", help="Skip CADD (faster)")
    parser.add_argument("--skip-civic", action="store_true", help="Skip CIViC")
    parser.add_argument("--from-tsv",   action="store_true",
                        help="Skip all API calls; load existing TSV and (re-)write VCF only")
    args = parser.parse_args()

    logger.info("=== Variant Annotation Pipeline ===")
    logger.info("Input:  %s", args.vcf)
    logger.info("Type:   %s", args.type)
    logger.info("Output: %s", args.output)

    # Checkpoint directory lives beside the output TSV
    ckpt_dir = Path(args.output).parent / ".annotation_cache"
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    ckpt = {i: str(ckpt_dir / f"layer{i}.json") for i in range(1, 6)}

    # ── Fast path: TSV already complete, just (re-)write VCF ──────────────────
    if args.from_tsv:
        if not Path(args.output).exists():
            logger.error("--from-tsv set but %s not found", args.output)
            sys.exit(1)
        logger.info("--from-tsv: loading %s, skipping all API calls", args.output)
        variants: list[dict] = []
        with open(args.output) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                variants.append(dict(row))
        logger.info("Loaded %d variants from TSV", len(variants))
        if args.vcf_output:
            write_annotated_vcf(variants, args.vcf, args.vcf_output)
        logger.info("=== Done (VCF-only mode) ===")
        return

    # ── Full annotation path with per-layer checkpointing ────────────────────

    # 1. Parse VCF (or resume from layer 1 checkpoint)
    if cached := load_checkpoint(ckpt[1]):
        variants = cached
    else:
        variants = parse_vcf(args.vcf, args.type)
        if not variants:
            logger.warning("No PASS variants found. Exiting.")
            sys.exit(0)

        # Layer 1: VEP
        logger.info("--- Layer 1: VEP ---")
        variants = VEPAnnotator().annotate(variants)
        logger.info("--- Normalizing HGVS ---")
        variants = [format_all(v) for v in variants]
        save_checkpoint(variants, ckpt[1])

    # Layer 2: gnomAD
    if cached := load_checkpoint(ckpt[2]):
        variants = cached
    else:
        logger.info("--- Layer 2: gnomAD ---")
        variants = GnomADAnnotator().annotate(variants)
        save_checkpoint(variants, ckpt[2])

    # Layer 3: ClinVar
    if cached := load_checkpoint(ckpt[3]):
        variants = cached
    else:
        logger.info("--- Layer 3: ClinVar ---")
        variants = ClinVarAnnotator().annotate(variants)
        save_checkpoint(variants, ckpt[3])

    # Layer 4: CADD
    if cached := load_checkpoint(ckpt[4]):
        variants = cached
    else:
        if not args.skip_cadd:
            logger.info("--- Layer 4: CADD ---")
            variants = CADDAnnotator().annotate(variants)
        else:
            for v in variants:
                v.update({"cadd_raw": ".", "cadd_phred": ".", "cadd_version": "."})
        save_checkpoint(variants, ckpt[4])

    # Layer 5: CIViC
    if cached := load_checkpoint(ckpt[5]):
        variants = cached
    else:
        if not args.skip_civic:
            logger.info("--- Layer 5: CIViC ---")
            variants = CIViCAnnotator().annotate(variants, variant_type=args.type)
        else:
            for v in variants:
                v.update({"civic_variant_id": ".", "civic_evidence_level": ".",
                          "civic_evidence_type": ".", "civic_clinical_significance": ".",
                          "civic_drug": ".", "civic_cancer_type": ".", "civic_evidence_count": "."})
        save_checkpoint(variants, ckpt[5])

    # Write outputs
    logger.info("--- Writing outputs ---")
    write_tsv(variants, args.output)
    if args.vcf_output:
        write_annotated_vcf(variants, args.vcf, args.vcf_output)

    # 9. Summary
    logger.info("=== Annotation complete ===")
    logger.info("Annotated %d variants", len(variants))

    high_impact = sum(1 for v in variants if v.get("impact") == "HIGH")
    pathogenic  = sum(1 for v in variants if "pathogenic" in v.get("clinvar_clnsig", "").lower())
    civic_hits  = sum(1 for v in variants if v.get("civic_evidence_level", ".") != ".")
    logger.info("  HIGH impact:       %d", high_impact)
    logger.info("  ClinVar Pathogenic:%d", pathogenic)
    logger.info("  CIViC evidence:    %d", civic_hits)

    if variants:
        logger.info("Sample annotated variant:")
        v = variants[0]
        logger.info("  %s", describe_variant(v))
        logger.info("  gnomAD AF: %s | CADD phred: %s | ClinVar: %s",
                    v.get("gnomad_af", "."),
                    v.get("cadd_phred", "."),
                    v.get("clinvar_clnsig", "."))


if __name__ == "__main__":
    main()
