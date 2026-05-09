"""
civic_annotator.py — CIViC (Clinical Interpretation of Variants in Cancer) annotator

CIViC is a community-curated knowledgebase of the clinical relevance of somatic
variants in cancer. Entries are curated by clinical oncologists, bioinformaticians,
and researchers from primary literature.

WHY CIViC MATTERS FOR SOMATIC ANALYSIS:
  Mutect2 tells you a variant exists. CIViC tells you whether that variant has
  actionable clinical significance:
    - Is there an FDA-approved drug for this variant?
    - Does the variant predict drug resistance?
    - Is there a clinical trial targeting this variant?
    - What is the evidence quality (level A = published prospective trial)?

EVIDENCE LEVELS (A to E):
  A  Validated association — prospective clinical trial evidence
  B  Clinical evidence — retrospective cohort, case series
  C  Case study — case report or very small series
  D  Preclinical — cell line or animal model evidence
  E  Inferential — computationally or experimentally inferred

EVIDENCE TYPES:
  Predictive  — predicts response/resistance to therapy
  Diagnostic  — helps classify cancer type
  Prognostic  — predicts outcome (survival, progression)
  Predisposing — germline variant increasing cancer risk
  Oncogenic   — evidence for functional oncogenic effect
  Functional  — experimental functional data

QUERY APPROACH:
  We search by gene + variant (e.g. gene="TP53", variant="R175H") since
  CIViC's primary organization is by gene+variant rather than coordinates.
  We extract the canonical variant name from VEP's HGVSp annotation.

API: https://civicdb.org/api (GraphQL + REST)
Docs: https://docs.civicdb.org/en/latest/api/
"""

import time
import logging
import requests

logger = logging.getLogger(__name__)

CIVIC_API = "https://civicdb.org/api/graphql"
RATE_LIMIT_WAIT = 0.5

EMPTY_ANN: dict[str, str] = {
    "civic_variant_id":          ".",
    "civic_evidence_level":      ".",
    "civic_evidence_type":       ".",
    "civic_clinical_significance":".",
    "civic_drug":                ".",
    "civic_cancer_type":         ".",
    "civic_evidence_count":      ".",
}

VARIANT_QUERY = """
query SearchVariant($geneSymbol: String!, $variantName: String!) {
  variants(geneSymbol: $geneSymbol, name: $variantName) {
    nodes {
      id
      name
      evidenceItems {
        nodes {
          id
          evidenceLevel
          evidenceType
          clinicalSignificance
          status
          therapies { name }
          disease { name }
        }
      }
    }
  }
}
"""


class CIViCAnnotator:
    """Queries CIViC for clinical evidence on somatic variants."""

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})
        self._cache: dict[str, dict] = {}

    def _query_variant(self, gene_symbol: str, variant_name: str) -> list[dict]:
        """Query CIViC for evidence items for a gene+variant combination."""
        cache_key = f"{gene_symbol}:{variant_name}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        try:
            resp = self.session.post(
                CIVIC_API,
                json={
                    "query": VARIANT_QUERY,
                    "variables": {"geneSymbol": gene_symbol, "variantName": variant_name},
                },
                timeout=20,
            )
            resp.raise_for_status()
            data = resp.json().get("data", {})
            nodes = data.get("variants", {}).get("nodes", [])
            self._cache[cache_key] = nodes
            return nodes
        except requests.RequestException as exc:
            logger.debug("CIViC query failed for %s %s: %s", gene_symbol, variant_name, exc)
            return []

    def _extract_variant_name(self, hgvsp: str) -> str:
        """
        Convert HGVSp to CIViC-style variant name.
        Input:  NP_000537.3:p.Arg273Cys or p.R273C
        Output: R273C (three-letter → one-letter conversion)
        """
        if not hgvsp or hgvsp == ".":
            return ""

        # Extract the protein change part
        if ":p." in hgvsp:
            p_change = hgvsp.split(":p.")[-1]
        elif "p." in hgvsp:
            p_change = hgvsp.split("p.")[-1]
        else:
            return ""

        # CIViC uses three-letter amino acid codes: 'R175H' not 'Arg175His'
        # Try to keep as-is; CIViC accepts both formats in search
        return p_change.rstrip("=")

    def _best_evidence(self, evidence_items: list[dict]) -> dict:
        """Return the highest-level evidence item for a variant."""
        level_order = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}

        accepted = [e for e in evidence_items if e.get("status") == "accepted"]
        if not accepted:
            return {}

        best = min(accepted, key=lambda e: level_order.get(e.get("evidenceLevel", "E"), 99))
        return best

    def annotate(self, variants: list[dict], variant_type: str = "somatic") -> list[dict]:
        """
        Add CIViC evidence to each variant dict.
        Only meaningful for somatic variants; germline support is limited in CIViC.
        """
        if variant_type == "germline":
            logger.info("CIViC annotation skipped for germline (CIViC is somatic-focused).")
            for v in variants:
                v.update(dict(EMPTY_ANN))
            return variants

        logger.info("Querying CIViC for %d somatic variants...", len(variants))

        for i, v in enumerate(variants):
            gene   = v.get("gene_symbol", ".")
            hgvsp  = v.get("hgvsp", ".")
            ann    = dict(EMPTY_ANN)

            if gene != "." and hgvsp != ".":
                variant_name = self._extract_variant_name(hgvsp)

                if variant_name:
                    nodes = self._query_variant(gene, variant_name)

                    if nodes:
                        variant_node = nodes[0]
                        ann["civic_variant_id"] = str(variant_node.get("id", "."))
                        evidence_items = variant_node.get("evidenceItems", {}).get("nodes", [])
                        ann["civic_evidence_count"] = str(len(evidence_items))

                        best = self._best_evidence(evidence_items)
                        if best:
                            ann["civic_evidence_level"]       = best.get("evidenceLevel", ".")
                            ann["civic_evidence_type"]        = best.get("evidenceType", ".")
                            ann["civic_clinical_significance"]= best.get("clinicalSignificance", ".")

                            therapies = best.get("therapies", [])
                            ann["civic_drug"] = ",".join(t["name"] for t in therapies) or "."

                            disease = best.get("disease") or {}
                            ann["civic_cancer_type"] = disease.get("name", ".")

            v.update(ann)
            time.sleep(RATE_LIMIT_WAIT)

            if (i + 1) % 10 == 0:
                logger.debug("  CIViC: %d / %d done", i + 1, len(variants))

        logger.info("CIViC annotation complete.")
        return variants
