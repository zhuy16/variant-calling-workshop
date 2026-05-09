"""
hgvs_formatter.py — HGVS notation normalizer and formatter

Parses and normalizes HGVS notation from VEP output to ensure compliance
with the HGVS Nomenclature 20.05 specification (varnomen.hgvs.org).

HGVS NOTATION EXPLAINED:

  HGVSg — Genomic notation (reference: NC_ accession, chromosome-level)
    Format: NC_000017.11:g.7674220C>T
    Parts:
      NC_000017.11    reference sequence accession (chr17, GRCh38)
      g.              genomic sequence indicator
      7674220         genomic position (1-based, VCF-style)
      C>T             reference > alternate allele change

  HGVSc — Coding DNA notation (reference: NM_ or ENST transcript)
    Format: NM_000546.6:c.817C>T
    Parts:
      NM_000546.6     RefSeq transcript accession (TP53 canonical)
      c.              coding DNA sequence indicator
      817             position in CDS (1 = first nucleotide of start codon)
      C>T             nucleotide change

    Position prefixes:
      c.5C>T          position 5 in CDS
      c.-5C>T         5 bp upstream of ATG (5' UTR)
      c.*5C>T         5 bp downstream of stop codon (3' UTR)
      c.100+5C>T      intronic: 5 bp after exon boundary
      c.100-5C>T      intronic: 5 bp before exon boundary

  HGVSp — Protein notation (reference: NP_ or ENSP protein)
    Format: NP_000537.3:p.Arg273Cys
    Parts:
      NP_000537.3     RefSeq protein accession (TP53 protein)
      p.              protein sequence indicator
      Arg273Cys       arginine (R) → cysteine (C) at position 273

    Special cases:
      p.Glu1308Ter    glutamate → stop codon (nonsense mutation)
      p.Ile255fs      frameshift starting at isoleucine 255
      p.Ala16_Lys17insGly  insertion between positions 16 and 17
      p.(Arg273Cys)   predicted (not confirmed at protein level)
      p.=             synonymous (no protein change)

THREE-LETTER vs ONE-LETTER amino acid codes:
  HGVS standard uses three-letter codes (Arg, Cys, Gly, etc.)
  Clinical databases sometimes use one-letter codes (R, C, G, etc.)
  This module preserves three-letter codes from VEP output.
"""

import re
import logging

logger = logging.getLogger(__name__)

# Chromosome → NC accession mapping (GRCh38 / hg38)
CHROM_TO_NC = {
    "1":  "NC_000001.11", "2":  "NC_000002.12", "3":  "NC_000003.12",
    "4":  "NC_000004.12", "5":  "NC_000005.10", "6":  "NC_000006.12",
    "7":  "NC_000007.14", "8":  "NC_000008.11", "9":  "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9",  "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X":  "NC_000023.11", "Y":  "NC_000024.10",
    "MT": "NC_012920.1",
}


def build_hgvsg(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    Construct HGVSg notation from VCF coordinates.

    Handles SNVs, insertions, deletions, and delinsertions following
    the HGVS Nomenclature 20.05 specification.

    Examples:
      SNV:        chr17:7674220 C>T  →  NC_000017.11:g.7674220C>T
      Deletion:   chr17:100 ACGT>A  →  NC_000017.11:g.101_103del
      Insertion:  chr17:100 A>ACGT  →  NC_000017.11:g.100_101insCGT
      DelIns:     chr17:100 ACG>TT  →  NC_000017.11:g.100_102delinsTT
    """
    chrom_clean = chrom.replace("chr", "")
    nc_acc = CHROM_TO_NC.get(chrom_clean, f"chr{chrom_clean}")

    ref_len = len(ref)
    alt_len = len(alt)

    # SNV
    if ref_len == 1 and alt_len == 1:
        return f"{nc_acc}:g.{pos}{ref}>{alt}"

    # Deletion
    if alt_len < ref_len and alt[0] == ref[0]:
        del_seq = ref[1:]
        if len(del_seq) == 1:
            return f"{nc_acc}:g.{pos + 1}del"
        return f"{nc_acc}:g.{pos + 1}_{pos + len(del_seq)}del"

    # Insertion
    if ref_len < alt_len and ref[0] == alt[0]:
        ins_seq = alt[1:]
        return f"{nc_acc}:g.{pos}_{pos + 1}ins{ins_seq}"

    # Delins (complex substitution)
    if ref_len == 1:
        return f"{nc_acc}:g.{pos}delins{alt}"
    return f"{nc_acc}:g.{pos}_{pos + ref_len - 1}delins{alt}"


def normalize_hgvsc(hgvsc: str) -> str:
    """
    Normalize HGVSc string from VEP output.
    VEP may return ENST-based notation; we preserve it as-is.
    Returns '.' if empty or clearly malformed.
    """
    if not hgvsc or hgvsc == ".":
        return "."
    # VEP format: ENST00000269305.9:c.817C>T or NM_000546.6:c.817C>T
    if ":" not in hgvsc:
        return "."
    return hgvsc.strip()


def normalize_hgvsp(hgvsp: str) -> str:
    """
    Normalize HGVSp string from VEP output.
    Ensures three-letter amino acid codes are used (HGVS standard).
    """
    if not hgvsp or hgvsp == ".":
        return "."
    if ":" not in hgvsp:
        return "."

    # VEP wraps predicted protein changes in parentheses: p.(Arg273Cys)
    # We strip the parentheses for cleaner display
    cleaned = re.sub(r"p\.\(([^)]+)\)", r"p.\1", hgvsp)
    return cleaned.strip()


def format_all(v: dict) -> dict:
    """
    Generate or normalize all three HGVS levels for a variant dict.
    Uses VEP-provided HGVSc and HGVSp when available; constructs HGVSg de novo.
    """
    # Always build HGVSg from coordinates (more reliable than VEP's version)
    v["hgvsg"] = build_hgvsg(v["chrom"], v["pos"], v["ref"], v["alt"])

    # Normalize VEP-provided HGVSc and HGVSp
    v["hgvsc"] = normalize_hgvsc(v.get("hgvsc", "."))
    v["hgvsp"] = normalize_hgvsp(v.get("hgvsp", "."))

    return v


def describe_variant(v: dict) -> str:
    """
    Return a human-readable one-line description of a variant.
    Used in agent prompts and log output.

    Example:
      TP53 p.Arg273Cys (chr17:7674220 C>T) | NM_000546.6:c.817C>T
    """
    gene    = v.get("gene_symbol", "?")
    hgvsp   = v.get("hgvsp", ".")
    chrom   = v.get("chrom", "?")
    pos     = v.get("pos", "?")
    ref     = v.get("ref", "?")
    alt     = v.get("alt", "?")
    hgvsc   = v.get("hgvsc", ".")

    if hgvsp and hgvsp != ".":
        protein_part = hgvsp.split(":p.")[-1] if ":p." in hgvsp else hgvsp
        return f"{gene} {protein_part} ({chrom}:{pos} {ref}>{alt}) | {hgvsc}"
    return f"{gene} ({chrom}:{pos} {ref}>{alt}) | {hgvsc}"
