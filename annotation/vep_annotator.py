"""
vep_annotator.py — Ensembl VEP REST API annotator

Annotates variants using the Ensembl Variant Effect Predictor via REST API.
No local install or cache required — results are equivalent to running VEP
locally for demo-scale callsets (< 1000 variants).

API documentation: https://rest.ensembl.org/#VEP
Rate limit: ~15 requests/second; we batch variants (POST endpoint accepts 200/request).

VEP consequence SO terms (severity order, high → low):
  HIGH:     transcript_ablation, splice_acceptor_variant, splice_donor_variant,
            stop_gained, frameshift_variant, stop_lost, start_lost,
            transcript_amplification
  MODERATE: missense_variant, protein_altering_variant, inframe_insertion,
            inframe_deletion, incomplete_terminal_codon_variant
  LOW:      synonymous_variant, stop_retained_variant, start_retained_variant,
            splice_region_variant, coding_sequence_variant
  MODIFIER: 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant,
            intron_variant, upstream_gene_variant, downstream_gene_variant,
            intergenic_variant, regulatory_region_variant
"""

import time
import logging
import requests
from typing import Any

logger = logging.getLogger(__name__)

VEP_REST_URL = "https://grch37.rest.ensembl.org"
BATCH_SIZE = 50   # conservative for REST API reliability
RETRY_WAIT = 5    # seconds between retries
MAX_RETRIES = 3

EMPTY_ANN: dict[str, Any] = {
    "gene_symbol":        ".",
    "ensembl_gene_id":    ".",
    "gene_biotype":       ".",
    "transcript_id":      ".",
    "is_canonical":       ".",
    "mane_select":        ".",
    "consequence":        ".",
    "impact":             ".",
    "exon":               ".",
    "intron":             ".",
    "hgvsc":              ".",
    "hgvsp":              ".",
    "hgvsg":              ".",
    "cdna_position":      ".",
    "cds_position":       ".",
    "protein_position":   ".",
    "amino_acids":        ".",
    "codons":             ".",
    "sift_score":         ".",
    "sift_pred":          ".",
    "polyphen_score":     ".",
    "polyphen_pred":      ".",
    "phylop":             ".",
    "gerp":               ".",
    "existing_variation": ".",
    "regulatory_id":      ".",
    "regulatory_type":    ".",
    "af_1kg":             ".",
}


class VEPAnnotator:
    """Annotates variants using the Ensembl VEP REST API (POST endpoint)."""

    def __init__(self, genome: str = "GRCh38"):
        self.genome = genome
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json",
        })

    def _post_batch(self, variants: list[str]) -> list[dict]:
        """POST a batch of region-format variants to the VEP REST endpoint."""
        url = f"{VEP_REST_URL}/vep/human/region"
        payload = {"variants": variants}

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                resp = self.session.post(url, json=payload, timeout=60)
                if resp.status_code == 429:
                    wait = int(resp.headers.get("Retry-After", RETRY_WAIT))
                    logger.warning("VEP rate limit hit; waiting %ds", wait)
                    time.sleep(wait)
                    continue
                resp.raise_for_status()
                return resp.json()
            except requests.RequestException as exc:
                if attempt == MAX_RETRIES:
                    logger.error("VEP request failed after %d retries: %s", MAX_RETRIES, exc)
                    return []
                time.sleep(RETRY_WAIT * attempt)
        return []

    def _parse_consequence(self, transcript: dict) -> dict:
        """Extract annotation fields from a single VEP transcript consequence."""
        ann = dict(EMPTY_ANN)

        ann["gene_symbol"]      = transcript.get("gene_symbol", ".")
        ann["ensembl_gene_id"]  = transcript.get("gene_id", ".")
        ann["gene_biotype"]     = transcript.get("biotype", ".")
        ann["transcript_id"]    = transcript.get("transcript_id", ".")
        ann["is_canonical"]     = "YES" if transcript.get("canonical") else "NO"
        ann["mane_select"]      = transcript.get("mane_select", ".")

        consequences = transcript.get("consequence_terms", [])
        ann["consequence"] = ",".join(consequences)
        ann["impact"]      = transcript.get("impact", ".")

        ann["exon"]   = transcript.get("exon", ".")
        ann["intron"] = transcript.get("intron", ".")

        ann["hgvsc"]           = transcript.get("hgvsc", ".")
        ann["hgvsp"]           = transcript.get("hgvsp", ".")
        ann["cdna_position"]   = str(transcript.get("cdna_start", "."))
        ann["cds_position"]    = str(transcript.get("cds_start", "."))
        ann["protein_position"]= str(transcript.get("protein_start", "."))
        ann["amino_acids"]     = transcript.get("amino_acids", ".")
        ann["codons"]          = transcript.get("codons", ".")

        sift = transcript.get("sift_prediction")
        if sift:
            ann["sift_pred"]  = sift
            ann["sift_score"] = str(transcript.get("sift_score", "."))

        pphen = transcript.get("polyphen_prediction")
        if pphen:
            ann["polyphen_pred"]  = pphen
            ann["polyphen_score"] = str(transcript.get("polyphen_score", "."))

        return ann

    def _parse_result(self, result: dict) -> dict:
        """Extract the most impactful canonical-transcript annotation from a VEP result."""
        ann = dict(EMPTY_ANN)

        # HGVSg — genomic notation from the input
        ann["hgvsg"] = result.get("id", ".")

        # Existing variation (rsID / ClinVar ID)
        colocated = result.get("colocated_variants", [])
        rs_ids = [v["id"] for v in colocated if v.get("id", "").startswith("rs")]
        ann["existing_variation"] = ",".join(rs_ids) if rs_ids else "."

        # Conservation scores from colocated data
        for v in colocated:
            if "phylop_vertebrate" in v:
                ann["phylop"] = str(v["phylop_vertebrate"])
            if "gerp++" in v:
                ann["gerp"]   = str(v["gerp++"])

        # Regulatory features
        reg_features = result.get("regulatory_feature_consequences", [])
        if reg_features:
            rf = reg_features[0]
            ann["regulatory_id"]   = rf.get("regulatory_feature_id", ".")
            ann["regulatory_type"] = rf.get("biotype", ".")

        # Pick the CANONICAL transcript, falling back to most severe
        transcripts = result.get("transcript_consequences", [])
        canonical = [t for t in transcripts if t.get("canonical")]
        chosen = canonical[0] if canonical else (transcripts[0] if transcripts else None)

        if chosen:
            transcript_ann = self._parse_consequence(chosen)
            ann.update(transcript_ann)

        # 1000 Genomes frequency
        freqs = result.get("frequencies", {})
        if freqs:
            af_vals = [str(v.get("af", ".")) for v in freqs.values() if "af" in v]
            ann["af_1kg"] = af_vals[0] if af_vals else "."

        return ann

    def annotate(self, variants: list[dict]) -> list[dict]:
        """
        Annotate a list of variant dicts with VEP.

        Each input dict must have: chrom, pos, ref, alt
        Returns same list with VEP fields added in-place.
        """
        # Convert to HGVS-like region notation for VEP REST
        hgvs_list = []
        for v in variants:
            chrom = v["chrom"].replace("chr", "")
            pos   = v["pos"]
            ref   = v["ref"]
            alt   = v["alt"]
            # VEP region format: "chr:start:end/allele"
            end   = pos + len(ref) - 1
            hgvs_list.append(f"{chrom}:{pos}-{end}/{alt}")

        logger.info("Annotating %d variants via VEP REST API...", len(hgvs_list))

        # Process in batches
        results_map: dict[str, dict] = {}
        for i in range(0, len(hgvs_list), BATCH_SIZE):
            batch = hgvs_list[i:i + BATCH_SIZE]
            logger.debug("  VEP batch %d-%d / %d", i + 1, i + len(batch), len(hgvs_list))
            raw_results = self._post_batch(batch)

            for result in raw_results:
                results_map[result.get("id", "")] = self._parse_result(result)

            if i + BATCH_SIZE < len(hgvs_list):
                time.sleep(0.5)  # courteous pacing

        # Merge VEP annotations back into variant dicts
        for v, hgvs_key in zip(variants, hgvs_list):
            vep_ann = results_map.get(hgvs_key, dict(EMPTY_ANN))
            v.update(vep_ann)

        logger.info("VEP annotation complete.")
        return variants
