#!/usr/bin/env python3
"""
variant_interpreter.py — Agentic variant interpretation layer

Takes the annotated TSV from the annotation pipeline and uses the Claude API
to generate clinical-grade natural language interpretations for each variant.

This mirrors a real tool a clinical bioinformatician would build to assist
with variant triage — converting structured annotation data into clear,
evidence-grounded clinical summaries.

WORKFLOW:
  Annotated TSV → Pre-classify by rule → Rank variants → Build structured
  prompt → Call Claude → Parse response → Write JSON output

PRE-CLASSIFICATION (rule-based, before LLM):
  Applied first to reduce unnecessary API calls. If the pre-classifier is
  confident (e.g., gnomAD AF > 5% → Benign), the LLM is skipped.

LLM CLASSIFICATION:
  For variants not pre-classified, sends structured prompt to Claude.
  Model: claude-opus-4-5 (temperature 0.2 for reproducibility)
  Demo mode: top N variants by impact + evidence score.

RANKING CRITERIA:
  Somatic:  impact score + CADD phred + CIViC evidence level + VAF
  Germline: impact score + CADD phred + ClinVar stars + gnomAD rarity

OUTPUT JSON FORMAT:
  {
    "metadata": { "pipeline_version", "model", "variant_type", "n_interpreted" },
    "variants": [
      {
        "variant_id":      "chr17:7674220:C:T",
        "gene":            "TP53",
        "hgvsp":           "p.Arg273Cys",
        "hgvsc":           "NM_000546.6:c.817C>T",
        "pre_classified":  false,
        "tier":            "DRIVER",
        "acmg_criteria":   ["PVS1"],          (germline only)
        "interpretation":  "This missense variant...",
        "evidence":        ["CADD phred 34", "COSMIC 12,000+ tumors", ...],
        "raw_annotation":  { ... all annotation fields ... }
      },
      ...
    ]
  }

Usage:
    export ANTHROPIC_API_KEY="sk-ant-..."
    python agent/variant_interpreter.py \\
        --tsv results/somatic/annotated.tsv \\
        --type somatic \\
        --top 10 \\
        --output results/somatic/interpreted.json
"""

import argparse
import csv
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import anthropic

# Load .env from project root if key not already in environment
if not os.environ.get("ANTHROPIC_API_KEY"):
    _env = Path(__file__).resolve().parent.parent / ".env"
    if _env.exists():
        for _line in _env.read_text().splitlines():
            if "=" in _line and not _line.startswith("#"):
                _k, _v = _line.split("=", 1)
                os.environ.setdefault(_k.strip(), _v.strip())

from prompts import (
    SYSTEM_PROMPT,
    build_somatic_prompt,
    build_germline_prompt,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

MODEL        = "claude-haiku-4-5-20251001"
MAX_TOKENS   = 600
TEMPERATURE  = 0.2
RATE_LIMIT_S = 1.0  # pause between API calls

IMPACT_SCORES = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
CIVIC_SCORES  = {"A": 4, "B": 3, "C": 2, "D": 1, "E": 0}


def load_tsv(tsv_path: str) -> list[dict]:
    """Load annotated variants TSV into list of dicts."""
    variants = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            variants.append(dict(row))
    logger.info("Loaded %d variants from %s", len(variants), tsv_path)
    return variants


def pre_classify(variant: dict, variant_type: str) -> str | None:
    """
    Rule-based pre-classification. Returns a tier string if confident,
    or None if the LLM should be consulted.

    Rules are deliberately conservative — the LLM handles ambiguous cases.
    """
    gnomad_af = _float(variant.get("gnomad_af"))
    clnsig    = (variant.get("clinvar_clnsig") or "").lower()
    impact    = (variant.get("impact") or "").upper()
    consequence = (variant.get("consequence") or "").lower()
    cadd      = _float(variant.get("cadd_phred"))

    if variant_type == "germline":
        # BA1: very common in gnomAD → Benign (stand-alone criterion)
        if gnomad_af is not None and gnomad_af > 0.05:
            return "BENIGN"
        # ClinVar Pathogenic with 3+ stars → high confidence
        stars = _int(variant.get("clinvar_review_stars"))
        if "pathogenic" in clnsig and "likely" not in clnsig and stars is not None and stars >= 3:
            return "PATHOGENIC"
        # Clearly benign ClinVar
        if "benign" in clnsig and "likely" not in clnsig and stars is not None and stars >= 2:
            return "BENIGN"
        # Synonymous with low CADD → BP7
        if "synonymous" in consequence and cadd is not None and cadd < 10:
            return "LIKELY_BENIGN"

    elif variant_type == "somatic":
        # Extremely common in gnomAD → likely germline contamination
        if gnomad_af is not None and gnomad_af > 0.01:
            return "ARTIFACT"
        # Modifier impact, no COSMIC evidence → likely passenger
        if impact == "MODIFIER" and variant.get("civic_evidence_level", ".") == ".":
            return "PASSENGER"

    return None


def score_variant(variant: dict, variant_type: str) -> float:
    """
    Compute a priority score for ranking variants before LLM interpretation.
    Higher = more likely to be clinically important = prioritized.
    """
    score = 0.0
    score += IMPACT_SCORES.get(variant.get("impact", ""), 0) * 10

    cadd = _float(variant.get("cadd_phred"))
    if cadd is not None:
        score += min(cadd, 50) * 0.5

    stars = _int(variant.get("clinvar_review_stars"))
    if stars is not None:
        score += stars * 5

    clnsig = (variant.get("clinvar_clnsig") or "").lower()
    if "pathogenic" in clnsig:
        score += 20
    if "benign" in clnsig:
        score -= 10

    if variant_type == "somatic":
        civic_level = variant.get("civic_evidence_level", ".")
        score += CIVIC_SCORES.get(civic_level, 0) * 8
        vaf = _float(variant.get("vaf"))
        if vaf is not None:
            score += vaf * 5  # higher VAF = clonal = more likely real driver

    if variant_type == "germline":
        gnomad_af = _float(variant.get("gnomad_af"))
        if gnomad_af is not None and gnomad_af < 0.0001:
            score += 15  # very rare → more likely pathogenic in rare disease context

    return score


def parse_llm_response(response_text: str, variant_type: str) -> dict:
    """Parse the structured LLM response into a dict."""
    result = {
        "tier": "UNCLASSIFIED",
        "interpretation": "",
        "acmg_criteria": [],
        "evidence": [],
        "raw_response": response_text,
    }

    lines = response_text.strip().split("\n")
    current_section = None
    evidence_lines: list[str] = []
    interp_lines:   list[str] = []

    for line in lines:
        line = line.strip()
        if line.startswith("INTERPRETATION:"):
            current_section = "interpretation"
            interp_lines.append(line[len("INTERPRETATION:"):].strip())
        elif line.startswith("TIER:"):
            tier_raw = line[len("TIER:"):].strip().upper()
            result["tier"] = tier_raw
            current_section = None
        elif line.startswith("CLASSIFICATION:"):
            cls_raw = line[len("CLASSIFICATION:"):].strip().upper()
            result["tier"] = cls_raw
            current_section = None
        elif line.startswith("ACMG_CRITERIA:"):
            criteria_str = line[len("ACMG_CRITERIA:"):].strip()
            result["acmg_criteria"] = [c.strip() for c in criteria_str.split(",") if c.strip()]
            current_section = None
        elif line.startswith("EVIDENCE:"):
            current_section = "evidence"
        elif current_section == "evidence" and line.startswith("-"):
            evidence_lines.append(line[1:].strip())
        elif current_section == "interpretation" and line:
            interp_lines.append(line)

    result["interpretation"] = " ".join(interp_lines).strip()
    result["evidence"]       = evidence_lines

    return result


def interpret_variant(variant: dict, client: anthropic.Anthropic,
                      variant_type: str, model: str = MODEL) -> dict:
    """Call Claude to interpret a single variant. Returns parsed response dict."""
    if variant_type == "somatic":
        prompt = build_somatic_prompt(variant)
    else:
        prompt = build_germline_prompt(variant)

    try:
        message = client.messages.create(
            model=model,
            max_tokens=MAX_TOKENS,
            temperature=TEMPERATURE,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": prompt}],
        )
        response_text = message.content[0].text
        return parse_llm_response(response_text, variant_type)
    except anthropic.APIError as exc:
        logger.error("Claude API error: %s", exc)
        return {
            "tier": "ERROR",
            "interpretation": f"API error: {exc}",
            "acmg_criteria": [],
            "evidence": [],
            "raw_response": str(exc),
        }


def _float(val: str | None) -> float | None:
    try:
        return float(val) if val and val != "." else None
    except (ValueError, TypeError):
        return None


def _int(val: str | None) -> int | None:
    try:
        return int(val) if val and val != "." else None
    except (ValueError, TypeError):
        return None


def main():
    parser = argparse.ArgumentParser(description="Agentic variant interpreter (Claude API)")
    parser.add_argument("--tsv",    required=True, help="Annotated variants TSV")
    parser.add_argument("--type",   required=True, choices=["somatic", "germline"])
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--top",    type=int, default=10,
                        help="Number of top variants to interpret via LLM (default: 10)")
    parser.add_argument("--all",    action="store_true",
                        help="Interpret all variants (overrides --top; careful with costs)")
    parser.add_argument("--model",   default=MODEL,
                        help="Claude model to use (default: claude-haiku-4-5-20251001)")
    args = parser.parse_args()

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        logger.error("ANTHROPIC_API_KEY environment variable not set.")
        logger.error("  export ANTHROPIC_API_KEY='sk-ant-...'")
        sys.exit(1)

    client = anthropic.Anthropic(api_key=api_key)

    # Load variants
    variants = load_tsv(args.tsv)
    if not variants:
        logger.error("No variants found in %s", args.tsv)
        sys.exit(1)

    # Pre-classify all variants by rule
    pre_classified_count = 0
    for v in variants:
        tier = pre_classify(v, args.type)
        if tier:
            v["_pre_tier"]       = tier
            v["_pre_classified"] = True
            pre_classified_count += 1
        else:
            v["_pre_tier"]       = None
            v["_pre_classified"] = False

    logger.info("Pre-classified %d / %d variants by rule", pre_classified_count, len(variants))

    # Rank all variants by priority score
    for v in variants:
        v["_score"] = score_variant(v, args.type)
    variants_sorted = sorted(variants, key=lambda v: v["_score"], reverse=True)

    # Determine which variants go to the LLM
    if args.all:
        llm_candidates = variants_sorted
    else:
        # Prioritize un-pre-classified, then fill to --top with pre-classified if needed
        unclassified = [v for v in variants_sorted if not v["_pre_classified"]]
        llm_candidates = unclassified[:args.top]

    logger.info("Will send %d variants to Claude (model=%s, temp=%.1f)",
                len(llm_candidates), MODEL, TEMPERATURE)

    # Estimate cost
    estimated_cost = len(llm_candidates) * 0.006  # rough estimate for claude-opus-4-5
    logger.info("Estimated API cost: ~$%.3f", estimated_cost)

    # Interpret each candidate
    results = []
    for i, v in enumerate(variants_sorted):
        variant_id = f"{v.get('chrom')}:{v.get('pos')}:{v.get('ref')}:{v.get('alt')}"

        if v["_pre_classified"]:
            result = {
                "variant_id":      variant_id,
                "gene":            v.get("gene_symbol", "."),
                "hgvsp":           v.get("hgvsp", "."),
                "hgvsc":           v.get("hgvsc", "."),
                "hgvsg":           v.get("hgvsg", "."),
                "consequence":     v.get("consequence", "."),
                "impact":          v.get("impact", "."),
                "pre_classified":  True,
                "tier":            v["_pre_tier"],
                "acmg_criteria":   [],
                "interpretation":  f"Pre-classified by rule: {v['_pre_tier']}",
                "evidence":        [],
                "priority_score":  round(v["_score"], 2),
                "raw_annotation":  {k: v[k] for k in v if not k.startswith("_")},
            }
        elif v in llm_candidates:
            logger.info("[%d/%d] Interpreting %s %s...",
                        i + 1, len(llm_candidates), v.get("gene_symbol", "?"),
                        v.get("hgvsp", v.get("consequence", "")))

            llm_result = interpret_variant(v, client, args.type, model=args.model)
            result = {
                "variant_id":     variant_id,
                "gene":           v.get("gene_symbol", "."),
                "hgvsp":          v.get("hgvsp", "."),
                "hgvsc":          v.get("hgvsc", "."),
                "hgvsg":          v.get("hgvsg", "."),
                "consequence":    v.get("consequence", "."),
                "impact":         v.get("impact", "."),
                "pre_classified": False,
                "tier":           llm_result["tier"],
                "acmg_criteria":  llm_result["acmg_criteria"],
                "interpretation": llm_result["interpretation"],
                "evidence":       llm_result["evidence"],
                "priority_score": round(v["_score"], 2),
                "raw_annotation": {k: v[k] for k in v if not k.startswith("_")},
            }
            time.sleep(RATE_LIMIT_S)
        else:
            # Variant not selected for LLM — include with score only
            result = {
                "variant_id":     variant_id,
                "gene":           v.get("gene_symbol", "."),
                "hgvsp":          v.get("hgvsp", "."),
                "hgvsc":          v.get("hgvsc", "."),
                "hgvsg":          v.get("hgvsg", "."),
                "consequence":    v.get("consequence", "."),
                "impact":         v.get("impact", "."),
                "pre_classified": False,
                "tier":           "NOT_INTERPRETED",
                "acmg_criteria":  [],
                "interpretation": "Below --top threshold; not sent to LLM.",
                "evidence":       [],
                "priority_score": round(v["_score"], 2),
                "raw_annotation": {k: v[k] for k in v if not k.startswith("_")},
            }

        results.append(result)

    # Write output JSON
    output = {
        "metadata": {
            "pipeline_version": "1.0.0",
            "model":            MODEL,
            "temperature":      TEMPERATURE,
            "variant_type":     args.type,
            "input_tsv":        args.tsv,
            "run_date":         datetime.now().isoformat(),
            "n_total":          len(variants),
            "n_pre_classified": pre_classified_count,
            "n_llm_interpreted":len(llm_candidates),
        },
        "variants": results,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(output, fh, indent=2)

    logger.info("Wrote %d interpreted variants to %s", len(results), args.output)

    # Print summary
    tier_counts: dict[str, int] = {}
    for r in results:
        tier = r["tier"]
        tier_counts[tier] = tier_counts.get(tier, 0) + 1

    logger.info("=== Tier Summary ===")
    for tier, count in sorted(tier_counts.items(), key=lambda x: -x[1]):
        logger.info("  %-25s %d", tier, count)

    # Print top 3 interpreted variants
    llm_results = [r for r in results if not r["pre_classified"] and r["tier"] != "NOT_INTERPRETED"]
    if llm_results:
        logger.info("\n=== Top Interpreted Variants ===")
        for r in llm_results[:3]:
            logger.info("  %s | %s | %s | %s",
                        r["gene"], r["hgvsp"], r["tier"],
                        r["interpretation"][:120] + "..." if len(r["interpretation"]) > 120
                        else r["interpretation"])


if __name__ == "__main__":
    main()
