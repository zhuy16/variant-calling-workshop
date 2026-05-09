#!/usr/bin/env bash
# =============================================================================
# somatic/03_filtering.sh
# Step 3 of 4 — Somatic Variant Filtering
#
# WHAT THIS STEP DOES:
#   1. Trains a read orientation artifact model from f1r2 counts
#   2. Applies FilterMutectCalls with all available evidence
#   3. Reports filter outcome statistics per variant
#
# THE FILTER CASCADE — what each Mutect2 filter catches:
#
#   clustered_events     — multiple somatic variants within a 35bp window.
#                         Real cancer variants can cluster (e.g. kataegis), but
#                         most clustered events are alignment artifacts near
#                         structural variant breakpoints.
#
#   contamination        — variant allele frequency is consistent with being
#                         a germline variant from a contaminating sample rather
#                         than a true somatic mutation.
#
#   duplicate_evidence   — only duplicate reads support this allele. Real
#                         somatic variants should have non-duplicate support.
#
#   fragment            — the forward and reverse reads in the same fragment
#                         disagree about the variant, which is inconsistent with
#                         a true variant.
#
#   germline             — the variant is likely germline based on its allele
#                         frequency in gnomAD and allele count in the normal.
#
#   haplotype            — variant is on a haplotype inconsistent with the
#                         assembled local sequence context.
#
#   low_allele_frac      — VAF below the minimum threshold (default 0.05).
#                         At this depth, very low VAF calls are unreliable.
#
#   map_qual             — reads supporting this allele have low mapping
#                         quality, suggesting they map to repetitive regions.
#
#   multiallelic         — site has more than 2 alleles — unusual for somatic
#                         SNPs; more often an alignment artifact at indels.
#
#   orientation          — the variant is enriched in reads with a specific
#                         orientation (F1R2 or F2R1) — the hallmark of FFPE
#                         oxidative damage (C>T or G>A in oxidized 8-oxoG).
#                         The LearnReadOrientationModel prior is applied here.
#
#   panel_of_normals     — variant appears in ≥2 normal samples in the PoN,
#                         indicating a recurrent sequencing/mapping artifact.
#
#   slippage             — strand slippage error near short tandem repeats.
#                         Polymerase slippage during PCR creates indels at
#                         microsatellite sites — a common false positive.
#
#   strand_bias          — statistically significant difference in allele
#                         support between forward and reverse reads.
#                         (Note: orientation is more specific for FFPE artifacts;
#                         strand_bias catches more generic biases.)
#
#   weak_evidence        — the log-odds score for somatic vs. artifact is
#                         below the threshold even after all other evidence.
#
# HOW FILTERS ARE APPLIED:
#   After FilterMutectCalls, each variant has a FILTER field:
#     PASS           — passed all filters; confident somatic variant
#     PASS + flags   — shouldn't happen; review if seen
#     filter_name    — failed one or more specific filters
#
#   Variants tagged PASS are the analysis-ready somatic callset.
#   You can still examine filtered variants to understand the biology
#   (e.g. many "contamination" calls in a sample may indicate real problems).
#
# INPUTS:  unfiltered.vcf.gz, f1r2.tar.gz, contamination.table, segments.table
# OUTPUTS: filtered.vcf.gz, filtering_stats.tsv
#
# RUN:     bash somatic/03_filtering.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

MUTECT_DIR="${ROOT}/results/somatic/mutect2"
OUTDIR="${ROOT}/results/somatic/filtered"
LOGDIR="${ROOT}/logs/somatic"
REF="${ROOT}/data/ref/Homo_sapiens_assembly38.chr17.fasta"
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] somatic/03 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/03_filtering.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)           VCF="$2";           shift 2 ;;
        --f1r2)          F1R2="$2";           shift 2 ;;
        --contamination) CONTAM_TABLE="$2";   shift 2 ;;
        --segments)      SEGMENTS="$2";       shift 2 ;;
        --outdir)        OUTDIR="$2";         shift 2 ;;
        *) shift ;;
    esac
done

VCF="${VCF:-${MUTECT_DIR}/unfiltered.vcf.gz}"
F1R2="${F1R2:-${MUTECT_DIR}/f1r2.tar.gz}"
CONTAM_TABLE="${CONTAM_TABLE:-${MUTECT_DIR}/contamination.table}"
SEGMENTS="${SEGMENTS:-${MUTECT_DIR}/segments.table}"

# ── Step 3A: Learn Read Orientation Model ─────────────────────────────────────
log "Learning read orientation artifact model from f1r2 counts..."
log "  This model captures FFPE oxidative damage: 8-oxoG lesions cause"
log "  C→T substitutions specifically in one strand orientation (F1R2 or F2R1)."
log "  The model assigns an orientation bias prior to each variant,   "
log "  which FilterMutectCalls then uses to downweight affected calls."

gatk --java-options "${MEM}" LearnReadOrientationModel \
    -I "${F1R2}" \
    -O "${OUTDIR}/artifact_priors.tar.gz" \
    2>&1 | tee -a "${LOGDIR}/03_filtering.log"

log "Orientation model trained. Priors saved to artifact_priors.tar.gz"

# ── Step 3B: FilterMutectCalls ────────────────────────────────────────────────
log ""
log "Applying FilterMutectCalls with full evidence set..."
log "  Contamination table: $(awk 'NR==2{printf \"%.4f\", $2}' "${CONTAM_TABLE}")"

gatk --java-options "${MEM}" FilterMutectCalls \
    -R "${REF}" \
    -V "${VCF}" \
    --contamination-table "${CONTAM_TABLE}" \
    --tumor-segmentation "${SEGMENTS}" \
    --ob-priors "${OUTDIR}/artifact_priors.tar.gz" \
    -O "${OUTDIR}/filtered.vcf.gz" \
    --filtering-stats "${OUTDIR}/filtering_stats.json" \
    2>&1 | tee -a "${LOGDIR}/03_filtering.log"

# ── Step 3C: Filter outcome summary ──────────────────────────────────────────
log ""
log "=== Filtering Summary ==="

TOTAL=$(bcftools view -H "${OUTDIR}/filtered.vcf.gz" | wc -l)
PASS=$(bcftools view -H -f PASS "${OUTDIR}/filtered.vcf.gz" | wc -l)
FILTERED=$((TOTAL - PASS))

log "  Total candidates:  ${TOTAL}"
log "  PASS (confident):  ${PASS}"
log "  Filtered out:      ${FILTERED}"
log "  Pass rate:         $(awk "BEGIN{printf \"%.1f%%\", ${PASS}/${TOTAL}*100}")"

log ""
log "  Filter breakdown (top reasons variants were removed):"
bcftools view "${OUTDIR}/filtered.vcf.gz" | \
    grep -v "^#" | \
    awk '{print $7}' | \
    grep -v PASS | \
    tr ";" "\n" | \
    sort | uniq -c | sort -rn | head -10 | \
    awk '{printf "  %-40s %d\n", $2, $1}' | \
    tee -a "${LOGDIR}/03_filtering.log"

log ""
log "  Saving filter breakdown to: ${OUTDIR}/filtering_stats.tsv"
echo -e "FILTER\tCOUNT" > "${OUTDIR}/filtering_stats.tsv"
bcftools view "${OUTDIR}/filtered.vcf.gz" | \
    grep -v "^#" | \
    awk '{print $7}' | \
    tr ";" "\n" | \
    sort | uniq -c | sort -rn | \
    awk '{print $2"\t"$1}' >> "${OUTDIR}/filtering_stats.tsv"

log ""
log "  INTERPRETING PASS RATE:"
log "  For chr17 tutorial data, expect ~50-200 PASS variants"
log "  Very high PASS counts (>1000) suggest filters may be too lenient"
log "  Very low PASS counts (<10) may indicate over-filtering or data issues"

log ""
log "=== Step 3 complete ==="
log "Outputs:"
log "  ${OUTDIR}/filtered.vcf.gz         — analysis-ready somatic callset"
log "  ${OUTDIR}/artifact_priors.tar.gz  — orientation bias model"
log "  ${OUTDIR}/filtering_stats.tsv     — per-filter counts"
log ""
log "Next step: bash somatic/04_annotation.sh"
