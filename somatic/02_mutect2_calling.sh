#!/usr/bin/env bash
# =============================================================================
# somatic/02_mutect2_calling.sh
# Step 2 of 4 — Mutect2 Tumor-Normal Calling + Contamination Estimation
#
# WHAT THIS STEP DOES:
#   1. Runs Mutect2 in tumor-normal mode to identify somatic variants
#   2. Collects pileup summaries at biallelic SNP sites
#   3. Estimates cross-sample contamination from the pileups
#
# HOW MUTECT2 WORKS (conceptually):
#   Mutect2 is not a simple threshold caller. It uses a local de novo
#   assembly approach:
#
#   1. ACTIVE REGION DETECTION: Identifies genomic windows with evidence
#      of variation (elevated mismatch rate vs. reference)
#
#   2. DE NOVO ASSEMBLY: Within each active region, assembles a De Bruijn
#      graph of observed k-mers (like a local genome assembler). This
#      generates candidate haplotypes that explain the reads.
#
#   3. LIKELIHOOD CALCULATION: For each candidate haplotype, calculates
#      the probability that the observed reads arose from:
#        - A somatic mutation (tumor allele + normal reference allele)
#        - A germline heterozygous variant (50% in both tumor and normal)
#        - A sequencing artifact
#
#   4. TUMOR-NORMAL COMPARISON: Uses the matched normal to:
#        - Filter germline variants (present at ~50% in both)
#        - Estimate baseline error rate for artifact calling
#        - Leverage the PoN for recurrent artifact sites
#
#   5. F1R2 READ ORIENTATION: Simultaneously counts reads by orientation
#      (F1R2 vs. F2R1) at each variant site. FFPE-damaged DNA shows
#      orientation bias (C→T mutations predominantly in one orientation).
#      This data feeds LearnReadOrientationModel in Step 3.
#
# CONTAMINATION ESTIMATION:
#   GetPileupSummaries measures allele counts at a set of common biallelic
#   SNPs from gnomAD. If the tumor has contamination from another sample,
#   you'll see unexpected alleles at sites where the tumor should be hom-ref.
#   CalculateContamination converts these counts into a contamination fraction.
#
#   Example output:
#     contamination = 0.02   (2% contamination — acceptable)
#     contamination = 0.08   (8% — high, needs careful filtering)
#     contamination = 0.25   (25% — likely sample swap or poor QC)
#
# INPUTS:  tumor.markdup.bam, normal.markdup.bam (from Step 1)
# OUTPUTS: unfiltered.vcf.gz, f1r2.tar.gz, contamination.table, segments.table
#
# RUN:     bash somatic/02_mutect2_calling.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

PREPDIR="${ROOT}/results/somatic/preprocessing"
OUTDIR="${ROOT}/results/somatic/mutect2"
LOGDIR="${ROOT}/logs/somatic"
REF="${ROOT}/data/ref/Homo_sapiens_assembly38.chr17.fasta"
GNOMAD="${ROOT}/data/resources/af-only-gnomad.chr17.vcf.gz"
PON="${ROOT}/data/resources/1000g_pon.chr17.vcf.gz"
INTERVALS="${ROOT}/data/resources/chr17_20M-60M.interval_list"
TUMOR_SAMPLE="HCC1143"
NORMAL_SAMPLE="HCC1143_NORMAL"
THREADS=4
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] somatic/02 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/02_mutect2.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --tumor-bam)     TUMOR_BAM="$2";    shift 2 ;;
        --normal-bam)    NORMAL_BAM="$2";   shift 2 ;;
        --tumor-sample)  TUMOR_SAMPLE="$2"; shift 2 ;;
        --normal-sample) NORMAL_SAMPLE="$2";shift 2 ;;
        --outdir)        OUTDIR="$2";       shift 2 ;;
        --threads)       THREADS="$2";      shift 2 ;;
        *) shift ;;
    esac
done

TUMOR_BAM="${TUMOR_BAM:-${PREPDIR}/tumor.markdup.bam}"
NORMAL_BAM="${NORMAL_BAM:-${PREPDIR}/normal.markdup.bam}"

# ── Step 2A: Mutect2 Tumor-Normal Calling ─────────────────────────────────────
log "Running Mutect2 (tumor-normal mode)..."
log "  Tumor:  ${TUMOR_SAMPLE}"
log "  Normal: ${NORMAL_SAMPLE}"
log "  Interval: chr17:7M-20M (TP53, BRCA1, ERBB2 region)"

gatk --java-options "${MEM}" Mutect2 \
    -R "${REF}" \
    -I "${TUMOR_BAM}" \
    -I "${NORMAL_BAM}" \
    -normal "${NORMAL_SAMPLE}" \
    --germline-resource "${GNOMAD}" \
    --panel-of-normals "${PON}" \
    -L "${INTERVALS}" \
    -O "${OUTDIR}/unfiltered.vcf.gz" \
    --f1r2-tar-gz "${OUTDIR}/f1r2.tar.gz" \
    --native-pair-hmm-threads "${THREADS}" \
    --max-mnp-distance 1 \
    2>&1 | tee -a "${LOGDIR}/02_mutect2.log"

log "Mutect2 done. Counting raw calls..."
NRAW=$(bcftools view -H "${OUTDIR}/unfiltered.vcf.gz" | wc -l)
log "  Raw somatic candidates: ${NRAW}"
log "  (This includes artifacts, germline leakage — filtering happens in Step 3)"

# ── Step 2B: GetPileupSummaries — Tumor ──────────────────────────────────────
log ""
log "Computing tumor pileup summaries for contamination estimation..."
log "  Measuring allele counts at common biallelic gnomAD SNPs (AF 0.05-0.95)"
log "  These sites let us detect unexpected alleles from sample contamination"

gatk --java-options "${MEM}" GetPileupSummaries \
    -I "${TUMOR_BAM}" \
    -V "${GNOMAD}" \
    -L "${INTERVALS}" \
    -O "${OUTDIR}/tumor_pileups.table" \
    2>&1 | tee -a "${LOGDIR}/02_mutect2.log"

# ── Step 2C: GetPileupSummaries — Normal ──────────────────────────────────────
log "Computing normal pileup summaries..."

gatk --java-options "${MEM}" GetPileupSummaries \
    -I "${NORMAL_BAM}" \
    -V "${GNOMAD}" \
    -L "${INTERVALS}" \
    -O "${OUTDIR}/normal_pileups.table" \
    2>&1 | tee -a "${LOGDIR}/02_mutect2.log"

# ── Step 2D: CalculateContamination ──────────────────────────────────────────
log "Estimating cross-sample contamination..."

gatk --java-options "${MEM}" CalculateContamination \
    -I "${OUTDIR}/tumor_pileups.table" \
    -matched "${OUTDIR}/normal_pileups.table" \
    -O "${OUTDIR}/contamination.table" \
    --tumor-segmentation "${OUTDIR}/segments.table" \
    2>&1 | tee -a "${LOGDIR}/02_mutect2.log"

log ""
log "=== Contamination Estimate ==="
CONTAM=$(awk 'NR==2{print $2}' "${OUTDIR}/contamination.table")
log "  Contamination fraction: ${CONTAM}"
log "  Interpretation:"
log "    < 0.02 = Low contamination (typical for well-QC'd samples)"
log "    0.02-0.05 = Moderate — FilterMutectCalls will compensate"
log "    > 0.05 = High — investigate sample handling; results less reliable"
log "    > 0.10 = Very high — possible sample swap; escalate to lab"

log ""
log "=== Step 2 complete ==="
log "Outputs:"
log "  ${OUTDIR}/unfiltered.vcf.gz    — raw Mutect2 calls (pre-filter)"
log "  ${OUTDIR}/f1r2.tar.gz          — read orientation counts (for Step 3)"
log "  ${OUTDIR}/contamination.table  — estimated contamination fraction"
log "  ${OUTDIR}/segments.table       — tumor copy number segments"
log ""
log "Next step: bash somatic/03_filtering.sh"
