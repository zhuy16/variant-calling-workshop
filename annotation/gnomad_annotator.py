"""
gnomad_annotator.py — gnomAD population allele frequency annotator

Queries the gnomAD GraphQL API for allele frequencies across 8 ancestry groups.
No local download required — the gnomAD v4 browser exposes a public GraphQL API.

gnomAD v4.1 dataset: 730,947 exomes + 76,215 genomes (GRCh38)
Ancestry groups:
  AFR — African / African American
  AMR — Latino / Admixed American
  ASJ — Ashkenazi Jewish
  EAS — East Asian
  FIN — Finnish
  NFE — Non-Finnish European
  SAS — South Asian
  MID — Middle Eastern

KEY CONCEPTS:
  AF (allele frequency): how common is the ALT allele across all sequenced individuals
  AC (allele count): raw count of ALT alleles observed
  AN (allele number): total alleles genotyped at this site (2 × samples with data)
  nhomalt: number of individuals homozygous for the ALT allele
  popmax_AF: highest AF across any single ancestry group (more conservative than global)
  faf95: filtering allele frequency at 95% confidence — the AF we can be 95% certain
         the variant doesn't exceed. Used in ACMG pathogenicity criteria (BA1: AF > 0.05).

INTERPRETING GNOMAD AF FOR CLINICAL INTERPRETATION:
  AF > 0.05 (5%)   → BA1: stand-alone Benign evidence (ACMG criterion)
                     Very unlikely to cause a severe recessive disease
  AF 0.001-0.05    → Common polymorphism; unlikely to be highly penetrant
  AF 0.0001-0.001  → Uncommon; may be relevant for dominant disease (>PM2 evidence)
  AF < 0.0001      → Rare; supports pathogenicity for dominant conditions (PM2)
  Not in gnomAD    → Very rare; strongest support for pathogenicity (PM2_supporting)

API: https://gnomad.broadinstitute.org/api (GraphQL)
"""

import time
import logging
import requests

logger = logging.getLogger(__name__)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
RATE_LIMIT_WAIT = 2.0  # seconds between requests (gnomAD enforces strict rate limits)

POPULATIONS = ["afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid"]

EMPTY_ANN: dict[str, str] = {
    "gnomad_af":          ".",
    "gnomad_ac":          ".",
    "gnomad_an":          ".",
    "gnomad_nhomalt":     ".",
    "gnomad_popmax_af":   ".",
    "gnomad_popmax_pop":  ".",
    "gnomad_faf95":       ".",
    "gnomad_afr_af":      ".",
    "gnomad_amr_af":      ".",
    "gnomad_asj_af":      ".",
    "gnomad_eas_af":      ".",
    "gnomad_fin_af":      ".",
    "gnomad_nfe_af":      ".",
    "gnomad_sas_af":      ".",
    "gnomad_mid_af":      ".",
    "gnomad_dataset":     ".",
}

QUERY = """
query VariantQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    genome {
      ac
      an
      af
      nhomalt
      faf95 { popmax popmax_population }
      populations {
        id
        ac
        an
        af
        nhomalt
      }
    }
    exome {
      ac
      an
      af
      nhomalt
      faf95 { popmax popmax_population }
      populations {
        id
        ac
        an
        af
        nhomalt
      }
    }
  }
}
"""


class GnomADAnnotator:
    """Queries gnomAD GraphQL API for population allele frequencies."""

    def __init__(self, dataset: str = "gnomad_r2_1"):
        self.dataset = dataset
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})

    def _query_variant(self, variant_id: str) -> dict:
        """
        Query gnomAD for a single variant.
        variant_id format: 'chrom-pos-ref-alt' (e.g. '17-7674220-C-T')
        """
        for attempt in range(4):
            try:
                resp = self.session.post(
                    GNOMAD_API,
                    json={"query": QUERY, "variables": {"variantId": variant_id, "dataset": self.dataset}},
                    timeout=30,
                )
                if resp.status_code == 429:
                    wait = 15 * (attempt + 1)
                    logger.warning("gnomAD 429 for %s, waiting %ds (attempt %d)", variant_id, wait, attempt + 1)
                    time.sleep(wait)
                    continue
                if not resp.ok:
                    logger.warning("gnomAD HTTP %d for %s", resp.status_code, variant_id)
                    return {}
                data = resp.json()
                return data.get("data", {}).get("variant") or {}
            except requests.RequestException as exc:
                logger.warning("gnomAD query failed for %s: %s", variant_id, exc)
                if attempt < 3:
                    time.sleep(5)
                    continue
                return {}
        logger.warning("gnomAD giving up on %s after 4 attempts", variant_id)
        return {}

    def _parse_frequencies(self, freq_data: dict | None) -> dict:
        """Parse genome or exome frequency block into annotation dict."""
        if not freq_data:
            return dict(EMPTY_ANN)

        ann = dict(EMPTY_ANN)
        ann["gnomad_af"]      = str(freq_data.get("af", ".") or ".")
        ann["gnomad_ac"]      = str(freq_data.get("ac", ".") or ".")
        ann["gnomad_an"]      = str(freq_data.get("an", ".") or ".")
        ann["gnomad_nhomalt"] = str(freq_data.get("nhomalt", ".") or ".")

        faf95 = freq_data.get("faf95") or {}
        ann["gnomad_faf95"]       = str(faf95.get("popmax", ".") or ".")
        ann["gnomad_popmax_af"]   = str(faf95.get("popmax", ".") or ".")
        ann["gnomad_popmax_pop"]  = str(faf95.get("popmax_population", ".") or ".")

        for pop in freq_data.get("populations", []):
            pop_id = pop.get("id", "").lower()
            if pop_id in POPULATIONS:
                ann[f"gnomad_{pop_id}_af"] = str(pop.get("af", ".") or ".")

        return ann

    def annotate(self, variants: list[dict]) -> list[dict]:
        """
        Add gnomAD population frequencies to each variant dict.
        Uses genome frequencies preferentially; falls back to exome.
        """
        logger.info("Querying gnomAD for %d variants...", len(variants))

        for i, v in enumerate(variants):
            chrom = v["chrom"].replace("chr", "")
            variant_id = f"{chrom}-{v['pos']}-{v['ref']}-{v['alt']}"

            result = self._query_variant(variant_id)

            # Prefer genome data; fall back to exome
            freq_data = result.get("genome") or result.get("exome")
            ann = self._parse_frequencies(freq_data)
            ann["gnomad_dataset"] = self.dataset

            v.update(ann)

            if (i + 1) % 100 == 0:
                hits = sum(1 for vv in variants[:i+1] if vv.get('gnomad_af', '.') != '.')
                logger.info("  gnomAD: %d/%d done, %d with AF data", i + 1, len(variants), hits)

            time.sleep(RATE_LIMIT_WAIT)

        logger.info("gnomAD annotation complete.")
        return variants
