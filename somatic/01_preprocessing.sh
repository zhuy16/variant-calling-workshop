#!/usr/bin/env bash
# =============================================================================
# somatic/01_preprocessing.sh
# Step 1 of 4 — Somatic Preprocessing
#
# WHAT THIS STEP DOES:
#   Prepares the tumor and normal BAMs for variant calling by removing
#   duplicate reads. This is the first mandatory step in every GATK workflow.
#
# WHY MARK DUPLICATES?
#   During library preparation, PCR amplification creates multiple copies of
#   the same DNA fragment. These copies are called "duplicates." In sequencing:
#
#   Real biology:          Read1 ──────── Read2      <- independent fragments
#   PCR duplicate:         Read1 ──────── Read2      <- copied from same template
#                          Read1 ──────── Read2      <- another copy
#
#   If duplicates are not removed:
#     - Apparent depth is artificially inflated (100× looks like 300×)
#     - Tumor VAF is unreliable: a variant at true 30% may appear at any %
#     - PCR errors introduced during amplification get counted multiple times
#       and can appear as false-positive somatic variants
#
#   GATK MarkDuplicates identifies duplicates by their 5' start positions:
#   pairs from the same original fragment will start at the same genomic
#   coordinate. It keeps ONE representative read and marks the rest as
#   duplicates in the FLAG field (0x400 bit). Downstream tools then ignore
#   flagged duplicates.
#
# OPTICAL DUPLICATES:
#   On patterned flowcells (NovaSeq, NextSeq 500), the array is physically
#   patterned — reads from adjacent array elements can be miscalled as
#   duplicates even without PCR. OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 is
#   the correct setting for patterned flowcells. Use 100 for HiSeq 2500.
#
# INPUTS:  tumor.bam, normal.bam (pre-aligned, from tutorial download)
# OUTPUTS: tumor.markdup.bam, normal.markdup.bam, per-sample metrics files
#
# RUN:     bash somatic/01_preprocessing.sh
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"
CFG="${ROOT}/config/pipeline_config.yml"

TUMOR_BAM="${ROOT}/data/bams/somatic/tumor.bam"
NORMAL_BAM="${ROOT}/data/bams/somatic/normal.bam"
OUTDIR="${ROOT}/results/somatic/preprocessing"
LOGDIR="${ROOT}/logs/somatic"
THREADS=4
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] somatic/01 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/01_preprocessing.log"; }

# ── Parse optional CLI args (used when called from Nextflow) ──────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)     TUMOR_BAM="$2";  shift 2 ;;
        --sample)  SAMPLE="$2";     shift 2 ;;
        --outdir)  OUTDIR="$2";     shift 2 ;;
        --threads) THREADS="$2";    shift 2 ;;
        *) shift ;;
    esac
done
SAMPLE="${SAMPLE:-tumor}"

# ── Step 1A: Mark Duplicates — Tumor ─────────────────────────────────────────
log "Marking duplicates in tumor BAM..."

gatk --java-options "${MEM}" MarkDuplicates \
    -I "${TUMOR_BAM}" \
    -O "${OUTDIR}/tumor.markdup.bam" \
    -M "${OUTDIR}/tumor.markdup.metrics" \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    2>&1 | tee -a "${LOGDIR}/01_preprocessing.log"

log "Tumor MarkDuplicates complete."

# ── Step 1B: Mark Duplicates — Normal ────────────────────────────────────────
log "Marking duplicates in normal BAM..."

gatk --java-options "${MEM}" MarkDuplicates \
    -I "${NORMAL_BAM}" \
    -O "${OUTDIR}/normal.markdup.bam" \
    -M "${OUTDIR}/normal.markdup.metrics" \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    2>&1 | tee -a "${LOGDIR}/01_preprocessing.log"

log "Normal MarkDuplicates complete."

# ── QC: Parse and display duplication metrics ─────────────────────────────────
log "=== Duplication Metrics Summary ==="
log "--- Tumor ---"
grep -A 2 "^LIBRARY" "${OUTDIR}/tumor.markdup.metrics" | tail -1 | \
    awk '{printf "  Estimated library size: %s\n  Percent duplication: %.2f%%\n", $9, $8*100}'
log "--- Normal ---"
grep -A 2 "^LIBRARY" "${OUTDIR}/normal.markdup.metrics" | tail -1 | \
    awk '{printf "  Estimated library size: %s\n  Percent duplication: %.2f%%\n", $9, $8*100}'

log ""
log "  INTERPRETING DUPLICATION RATE:"
log "  < 10%  = Low duplication (excellent library, high complexity)"
log "  10-25% = Acceptable for WGS/WES"
log "  25-50% = High — consider increasing input DNA or reducing PCR cycles"
log "  > 50%  = Very high — library likely over-amplified; depth is misleading"

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: BWA-MEM Alignment                                             ║
# ║                                                                          ║
# ║  WHAT IT DOES: Maps raw FASTQ reads to the reference genome.             ║
# ║  BWA-MEM2 aligns 150bp paired-end reads using a BWT-FM index of         ║
# ║  the reference, outputting an aligned BAM sorted by coordinate.          ║
# ║                                                                          ║
# ║  WHY SKIPPED: The GATK tutorial BAMs are pre-aligned. Starting from      ║
# ║  BAM is correct for this tutorial. In production, you start from         ║
# ║  FASTQ (output of the sequencer's BCL2FASTQ demultiplexer).              ║
# ║                                                                          ║
# ║  RUNTIME: ~4 hours per 30× WGS sample on 32 threads.                    ║
# ║  MEM: ~8 GB RAM for 3 GB hg38 index.                                    ║
# ║                                                                          ║
# ║  EXPECTED OUTPUT: aligned.bam, ~45 GB for 30× WGS                       ║
# ║  samtools flagstat should show:                                          ║
# ║    99.5%+ reads mapped                                                   ║
# ║    95%+   properly paired                                                ║
# ║    < 1%   singletons                                                     ║
# ║                                                                          ║
# ║  TO ENABLE: Set FASTQ_R1 and FASTQ_R2, uncomment below.                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# FASTQ_R1="data/fastq/tumor_R1.fastq.gz"
# FASTQ_R2="data/fastq/tumor_R2.fastq.gz"
# READ_GROUP="@RG\tID:tumor\tSM:HCC1143\tPL:ILLUMINA\tLB:lib1\tPU:flowcell1"
#
# bwa mem -t ${THREADS} -R "${READ_GROUP}" \
#     data/ref/Homo_sapiens_assembly38.chr17.fasta \
#     ${FASTQ_R1} ${FASTQ_R2} \
#     | samtools sort -@ ${THREADS} -o aligned_tumor.bam
# samtools index aligned_tumor.bam

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: BQSR (Base Quality Score Recalibration)                        ║
# ║                                                                          ║
# ║  WHAT IT DOES: The sequencer reports a quality score (Phred Q) for each  ║
# ║  base call — e.g. Q30 means 1-in-1000 error probability. But these       ║
# ║  scores have systematic biases: the sequencer's reported Q30 may         ║
# ║  actually be Q28 for certain base contexts (e.g. GGT triplets) due to    ║
# ║  known chemistry effects. BQSR models these biases using known variant    ║
# ║  sites and corrects the reported quality scores across the BAM.           ║
# ║                                                                          ║
# ║  COVARIATES MODELED:                                                     ║
# ║    1. Read group (sequencing run)                                        ║
# ║    2. Reported quality score                                             ║
# ║    3. Machine cycle (position in read, e.g. end of read is lower Q)     ║
# ║    4. Dinucleotide context (e.g. GA context has inflated error rate)     ║
# ║                                                                          ║
# ║  EFFECT ON VARIANT CALLING: Miscalibrated base qualities can cause       ║
# ║  HaplotypeCaller/Mutect2 to over-call or under-call variants, because    ║
# ║  the variant likelihood model uses base quality directly.                ║
# ║                                                                          ║
# ║  WHY SKIPPED: Requires 16–32 GB RAM on full chr-scale BAMs.              ║
# ║  On modern NovaSeq (patterned flowcell), base quality is well-           ║
# ║  calibrated out of the box — BQSR benefit is 0.1–0.5% sensitivity on    ║
# ║  these platforms. The GATK team acknowledges this in their docs.         ║
# ║                                                                          ║
# ║  EXPECTED OUTPUT: recal.table shows per-covariate error rates.           ║
# ║  AnalyzeCovariates plot should show Q-score distribution shift:          ║
# ║    Before: Q30 actual error rate = 0.12% (slightly miscalibrated)       ║
# ║    After:  Q30 actual error rate = 0.10% (matches Phred expectation)    ║
# ║                                                                          ║
# ║  TO ENABLE: uncomment below; set MEM to at least 16g                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# for SAMPLE_BAM in "${OUTDIR}/tumor.markdup.bam" "${OUTDIR}/normal.markdup.bam"; do
#     BASE=$(basename "${SAMPLE_BAM}" .markdup.bam)
#     gatk --java-options "-Xmx16g" BaseRecalibrator \
#         -I "${SAMPLE_BAM}" \
#         -R data/ref/Homo_sapiens_assembly38.chr17.fasta \
#         --known-sites data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf.gz \
#         --known-sites data/resources/Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf.gz \
#         -O "${OUTDIR}/${BASE}.recal.table"
#
#     gatk --java-options "-Xmx16g" ApplyBQSR \
#         -I "${SAMPLE_BAM}" \
#         -R data/ref/Homo_sapiens_assembly38.chr17.fasta \
#         --bqsr-recal-file "${OUTDIR}/${BASE}.recal.table" \
#         -O "${OUTDIR}/${BASE}.bqsr.bam"
# done

log ""
log "=== Step 1 complete ==="
log "Outputs:"
log "  ${OUTDIR}/tumor.markdup.bam"
log "  ${OUTDIR}/normal.markdup.bam"
log "  ${OUTDIR}/tumor.markdup.metrics"
log "  ${OUTDIR}/normal.markdup.metrics"
log ""
log "Next step: bash somatic/02_mutect2_calling.sh"
