"""
cadd_annotator.py — CADD (Combined Annotation Dependent Depletion) score annotator

CADD is a machine learning framework that integrates 60+ annotation features
(conservation, regulatory, transcript impact, epigenomic data, etc.) into a
single deleteriousness score for any SNV or indel in the human genome.

HOW CADD WORKS:
  CADD trains a support vector machine to distinguish:
    DERIVED alleles (evolved substitutions that were tolerated by natural selection)
    vs.
    SIMULATED alleles (random mutations from a neutral model)
  Variants that look more like simulated mutations score higher — they are more
  likely to be deleterious because evolution never accepted them.

THE TWO CADD SCORES:
  CADD raw:   raw SVM score (negative = evolved/tolerated; positive = simulated/deleterious)
  CADD phred: scaled to a Phred-like format: CADD_phred = -10 * log10(rank / total)
              This makes interpretation intuitive:
                CADD phred 10 = top 10% most deleterious in the genome
                CADD phred 20 = top 1%
                CADD phred 30 = top 0.1%
                CADD phred 40 = top 0.01% (extremely rare to see in healthy individuals)

INTERPRETING CADD IN CLINICAL CONTEXT:
  For germline variants (ACMG framework):
    phred ≥ 25 supports PP3 (in silico evidence of pathogenicity)
    phred < 15 supports BP4 (in silico evidence of benign impact)
  Note: CADD alone is NOT sufficient to classify — always combine with other evidence.

  For somatic variants:
    High CADD phred scores in driver genes (TP53, KRAS, etc.) reinforce pathogenicity
    Low CADD scores in non-canonical positions suggest passenger mutations

CADD API:
  Free, no registration: https://cadd.gs.washington.edu/api
  Correct endpoint (GET, returns JSON):
    /api/v1.0/<version>/<chrom>:<pos>-<pos>
  Example:
    curl https://cadd.gs.washington.edu/api/v1.0/GRCh37-v1.6/20:9999996-9999996
"""

import time
import logging
import requests

logger = logging.getLogger(__name__)

CADD_API = "https://cadd.gs.washington.edu/api/v1.0"
RATE_LIMIT_WAIT = 0.25  # GET requests are lightweight; be courteous

EMPTY_ANN: dict[str, str] = {
    "cadd_raw":   ".",
    "cadd_phred": ".",
    "cadd_version": ".",
}


class CADDAnnotator:
    """Retrieves CADD scores via the CADD web API (GRCh37, v1.6).

    Uses GET /api/v1.0/{version}/{chrom}:{pos}-{pos} which returns a JSON
    array of all three possible SNVs at that position.  We match by ref+alt.
    """

    def __init__(self, genome: str = "GRCh37-v1.6"):
        self.genome = genome
        self.session = requests.Session()

    def _lookup(self, chrom: str, pos: str) -> list[dict]:
        """
        GET all SNVs at a single position.
        Returns a list of dicts with keys: Chrom, Pos, Ref, Alt, RawScore, PHRED
        """
        url = f"{CADD_API}/{self.genome}/{chrom}:{pos}-{pos}"
        try:
            resp = self.session.get(url, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            # Response is [[header...], [row...], ...] OR [] if no pre-computed score
            if not data or not isinstance(data[0], list):
                return []
            header = [h.lower() for h in data[0]]
            return [dict(zip(header, row)) for row in data[1:]]
        except Exception as exc:
            logger.warning("CADD lookup failed %s:%s — %s", chrom, pos, exc)
            return []

    def annotate(self, variants: list[dict]) -> list[dict]:
        """Add CADD raw and phred scores to each variant dict."""
        logger.info("Fetching CADD scores for %d variants...", len(variants))
        hits = 0

        for idx, v in enumerate(variants):
            chrom = str(v["chrom"]).replace("chr", "")
            pos   = str(v["pos"])
            ref   = str(v["ref"]).upper()
            alt   = str(v["alt"]).upper()

            rows = self._lookup(chrom, pos)
            score = next(
                (r for r in rows
                 if str(r.get("ref", "")).upper() == ref
                 and str(r.get("alt", "")).upper() == alt),
                None,
            )

            if score:
                v["cadd_raw"]   = str(score.get("rawscore", "."))
                v["cadd_phred"] = str(score.get("phred",    "."))
                hits += 1
            else:
                v["cadd_raw"]   = "."
                v["cadd_phred"] = "."
            v["cadd_version"] = self.genome

            time.sleep(RATE_LIMIT_WAIT)

            if (idx + 1) % 100 == 0:
                logger.info("  CADD: %d/%d done, %d scored", idx + 1, len(variants), hits)

        logger.info("CADD annotation complete. %d/%d variants scored.", hits, len(variants))
        return variants
