"""
clinvar_annotator.py — NCBI ClinVar annotation via eUtils API

Retrieves clinical significance and disease associations from ClinVar using
NCBI's eUtils API (no registration, no API key required for low-volume queries).

ClinVar (clinicalgenome.org): the primary public resource for variant-disease
interpretations. Submissions come from clinical labs, research groups, and
expert panels. The star rating reflects the evidence quality:

  ★★★★ (4 stars): Reviewed by expert panel (highest reliability)
  ★★★  (3 stars): Reviewed by expert panel (provisional)
  ★★   (2 stars): Criteria provided, multiple submitters, no conflicts
  ★    (1 star):  Criteria provided, single submitter
  ☆    (0 stars): No assertion criteria provided (lowest reliability)

CLINICAL SIGNIFICANCE CATEGORIES:
  Pathogenic (P)          — variant causes disease
  Likely pathogenic (LP)  — variant probably causes disease (≥90% probability)
  Uncertain significance (VUS) — insufficient or conflicting evidence
  Likely benign (LB)      — variant probably does not cause disease
  Benign (B)              — variant does not cause disease

ACMG CLASSIFICATION NOTE:
  ClinVar stores the classifications from labs, but classification itself
  follows the ACMG/AMP 2015 guidelines. Key criteria:
    PVS1 — Loss-of-function in a gene where LOF is a known mechanism
    PS1   — Same amino acid change as an established Pathogenic variant
    PM2   — Absent from gnomAD (or very rare) in a recessive gene
    BA1   — Allele frequency > 5% in gnomAD → stand-alone Benign evidence

API: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
"""

import time
import logging
import requests

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RATE_LIMIT_WAIT = 0.34  # NCBI allows 3 requests/sec without API key

EMPTY_ANN: dict[str, str] = {
    "clinvar_clnsig":       ".",
    "clinvar_clndn":        ".",
    "clinvar_review_stars": ".",
    "clinvar_review_status":".",
    "clinvar_allele_id":    ".",
    "clinvar_variation_id": ".",
    "clinvar_last_evaluated":".",
    "clinvar_submission_count": ".",
}

STAR_MAP = {
    "no_assertion_provided":                   "0",
    "no_assertion_criteria_provided":          "0",
    "criteria_provided,_single_submitter":     "1",
    "criteria_provided,_conflicting_interpretations": "1",
    "criteria_provided,_multiple_submitters,_no_conflicts": "2",
    "reviewed_by_expert_panel":                "3",
    "practice_guideline":                      "4",
}


class ClinVarAnnotator:
    """Retrieves ClinVar classifications using NCBI eUtils (esearch + esummary)."""

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    def _esearch(self, chrom: str, pos: int, ref: str, alt: str) -> str | None:
        """Search ClinVar by genomic coordinates. Returns ClinVar variation ID or None."""
        # ClinVar esearch accepts HGVS-like queries
        query = f"{chrom}[chr] AND {pos}[chrpos37] AND {ref}>{alt}[varnt]"
        try:
            resp = self.session.get(
                f"{EUTILS_BASE}/esearch.fcgi",
                params={"db": "clinvar", "term": query, "retmode": "json", "retmax": 1},
                timeout=15,
            )
            resp.raise_for_status()
            ids = resp.json().get("esearchresult", {}).get("idlist", [])
            return ids[0] if ids else None
        except requests.RequestException as exc:
            logger.debug("ClinVar esearch failed for %s:%d: %s", chrom, pos, exc)
            return None

    def _esummary(self, variation_id: str) -> dict:
        """Fetch ClinVar summary for a variation ID."""
        try:
            resp = self.session.get(
                f"{EUTILS_BASE}/esummary.fcgi",
                params={"db": "clinvar", "id": variation_id, "retmode": "json"},
                timeout=15,
            )
            resp.raise_for_status()
            result = resp.json().get("result", {})
            return result.get(variation_id, {})
        except requests.RequestException as exc:
            logger.debug("ClinVar esummary failed for %s: %s", variation_id, exc)
            return {}

    def _parse_summary(self, summary: dict, variation_id: str) -> dict:
        """Parse a ClinVar esummary result into annotation fields."""
        ann = dict(EMPTY_ANN)
        ann["clinvar_variation_id"] = variation_id

        # Clinical significance
        clinsig = summary.get("clinical_significance", {})
        ann["clinvar_clnsig"]          = clinsig.get("description", ".")
        ann["clinvar_last_evaluated"]  = clinsig.get("last_evaluated", ".")
        ann["clinvar_submission_count"]= str(clinsig.get("submission_count", "."))

        # Review status → star rating
        review_status = clinsig.get("review_status", ".")
        ann["clinvar_review_status"] = review_status
        normalized = review_status.lower().replace(" ", "_")
        ann["clinvar_review_stars"]  = STAR_MAP.get(normalized, "0")

        # Disease name
        trait_set = summary.get("trait_set", [])
        if trait_set:
            trait_names = [t.get("trait_name", "") for t in trait_set if t.get("trait_name")]
            ann["clinvar_clndn"] = "|".join(trait_names) if trait_names else "."

        # Allele ID
        allele_id = summary.get("allele_id")
        if allele_id:
            ann["clinvar_allele_id"] = str(allele_id)

        return ann

    def annotate(self, variants: list[dict]) -> list[dict]:
        """Add ClinVar annotations to each variant dict."""
        logger.info("Querying ClinVar for %d variants...", len(variants))

        for i, v in enumerate(variants):
            chrom = v["chrom"].replace("chr", "")
            variation_id = self._esearch(chrom, v["pos"], v["ref"], v["alt"])

            if variation_id:
                summary = self._esummary(variation_id)
                ann = self._parse_summary(summary, variation_id)
            else:
                ann = dict(EMPTY_ANN)

            v.update(ann)

            if (i + 1) % 10 == 0:
                logger.debug("  ClinVar: %d / %d done", i + 1, len(variants))

            time.sleep(RATE_LIMIT_WAIT)

        logger.info("ClinVar annotation complete.")
        return variants
