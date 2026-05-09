#!/usr/bin/env bash
# =============================================================================
# germline/01_preprocessing.sh
# Step 1 of 4 — Germline Preprocessing (trio: mother, father, son)
#
# WHAT THIS STEP DOES:
#   Marks duplicate reads in all three trio member BAMs, preparing them for
#   per-sample HaplotypeCaller GVCF generation in Step 2.
#
# WHY A TRIO?
#   Joint calling across a mother-father-son trio enables:
#
#   1. DE NOVO DETECTION: A variant in the son absent from both parents is
#      a de novo mutation. These are important in rare pediatric disease:
#      ~50% of severe neurodevelopmental disorders are caused by de novos.
#      Detection requires all three samples to be called jointly.
#
#   2. COMPOUND HET PHASING: If a gene has two heterozygous variants (one
#      on each allele), the parental genotypes reveal which are cis (same
#      allele — usually benign) vs. trans (different alleles — potentially
#      pathogenic in autosomal recessive disease).
#
#   3. MENDELIAN CONSISTENCY: A child can only carry alleles inherited from
#      its parents. Variants that violate Mendelian inheritance (neither
#      parent carries the allele) are flagged as de novo candidates or
#      potential sequencing errors — a built-in quality metric.
#
#   4. STATISTICAL POWER: Even at moderate depth, joint genotyping across
#      trio members gives more accurate genotype calls than single-sample
#      calling, especially for low-coverage or repetitive regions.
#
# THE GATK TUTORIAL TRIO:
#   NA12878 (mother) — the Genome in a Bottle (GIAB) benchmark sample.
#     Over 4 million variants have been truth-validated in this genome.
#   NA12891 (father) — HapMap Caucasian reference individual
#   NA12892 (son)    — child of NA12878 + NA12891
#
# INPUTS:  mother.bam, father.bam, son.bam (chr20 subset, ~20 MB each)
# OUTPUTS: *.markdup.bam for each trio member
#
# RUN:     bash germline/01_preprocessing.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

BAM_DIR="${ROOT}/data/bams/germline"
OUTDIR="${ROOT}/results/germline/preprocessing"
LOGDIR="${ROOT}/logs/germline"
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] germline/01 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/01_preprocessing.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)     INPUT_BAM="$2"; shift 2 ;;
        --sample)  SAMPLE="$2";    shift 2 ;;
        --outdir)  OUTDIR="$2";    shift 2 ;;
        --threads) THREADS="$2";   shift 2 ;;
        *) shift ;;
    esac
done

# ── If called directly (not from Nextflow), process all three trio members ────
if [[ -z "${INPUT_BAM:-}" ]]; then
    for SAMPLE in mother father son; do
        INPUT_BAM="${BAM_DIR}/${SAMPLE}.bam"
        if [[ ! -f "${INPUT_BAM}" ]]; then
            log "ERROR: ${INPUT_BAM} not found. Run: bash data/download_tutorial_data.sh"
            exit 1
        fi

        log "Processing ${SAMPLE}..."
        gatk --java-options "${MEM}" MarkDuplicates \
            -I "${INPUT_BAM}" \
            -O "${OUTDIR}/${SAMPLE}.markdup.bam" \
            -M "${OUTDIR}/${SAMPLE}.markdup.metrics" \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            2>&1 | tee -a "${LOGDIR}/01_preprocessing.log"

        DUP_RATE=$(grep -A 2 "^LIBRARY" "${OUTDIR}/${SAMPLE}.markdup.metrics" | \
            tail -1 | awk '{printf "%.1f%%", $8*100}')
        log "  ${SAMPLE}: duplication rate = ${DUP_RATE}"
    done
else
    # Called from Nextflow for a single sample
    SAMPLE="${SAMPLE:-unknown}"
    log "Processing single sample: ${SAMPLE}"
    gatk --java-options "${MEM}" MarkDuplicates \
        -I "${INPUT_BAM}" \
        -O "${OUTDIR}/${SAMPLE}.markdup.bam" \
        -M "${OUTDIR}/${SAMPLE}.markdup.metrics" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        2>&1 | tee -a "${LOGDIR}/01_preprocessing.log"
fi

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: Read QC (FastQC / MultiQC)                                    ║
# ║                                                                          ║
# ║  WHAT IT DOES: Generates a comprehensive HTML quality report for each    ║
# ║  FASTQ file, covering: per-base quality scores, sequence length          ║
# ║  distribution, GC content bias, adapter contamination, and              ║
# ║  over-represented sequences.                                             ║
# ║                                                                          ║
# ║  WHY SKIPPED: Tutorial BAMs are pre-QC'd synthetic data. In production,  ║
# ║  run FastQC on raw FASTQs before ANY other step — this is the gate that  ║
# ║  catches failed lanes, wrong samples, and adapter contamination.         ║
# ║                                                                          ║
# ║  RED FLAGS to look for in FastQC:                                        ║
# ║    - Per-base quality drops below Q20 toward read ends (normal for       ║
# ║      most platforms, trim if severe)                                     ║
# ║    - GC content bimodal distribution (possible contamination)            ║
# ║    - Adapter content > 5% (trim with fastp or Trimmomatic)               ║
# ║    - Per-sequence quality peak < Q28 (bad run)                           ║
# ║                                                                          ║
# ║  EXPECTED: Most bases Q30+, GC ~50%, adapters < 1%, duplication 10-25%  ║
# ║                                                                          ║
# ║  TO ENABLE: pip install multiqc; uncomment below                         ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# for FASTQ in data/fastq/*.fastq.gz; do
#     fastqc "${FASTQ}" -o results/germline/fastqc/ -t 4
# done
# multiqc results/germline/fastqc/ -o results/germline/multiqc/

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: Adapter Trimming (fastp / Trimmomatic)                         ║
# ║                                                                          ║
# ║  WHAT IT DOES: Removes Illumina adapter sequences and low-quality bases  ║
# ║  from the ends of reads before alignment.                                ║
# ║                                                                          ║
# ║  WHY ADAPTERS MATTER: If adapter sequence is present in reads, BWA will  ║
# ║  either fail to align them (losing coverage) or mis-align them with      ║
# ║  soft-clipping (distorting variant calls near read ends).                ║
# ║                                                                          ║
# ║  WHY SKIPPED: Tutorial data is pre-trimmed. fastp is preferred for modern║
# ║  data — it handles adapter auto-detection and quality trimming in one    ║
# ║  pass at high speed.                                                     ║
# ║                                                                          ║
# ║  TO ENABLE: conda install -c bioconda fastp; uncomment below             ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# fastp \
#     -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
#     -o sample_R1_trimmed.fastq.gz -O sample_R2_trimmed.fastq.gz \
#     --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#     --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     --cut_right --cut_right_mean_quality 20 \
#     --json results/fastqc/sample_fastp.json \
#     --html results/fastqc/sample_fastp.html \
#     --thread 4

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: BQSR (Base Quality Score Recalibration)                        ║
# ║                                                                          ║
# ║  Same as somatic module — see somatic/01_preprocessing.sh for full       ║
# ║  explanation. For germline calling, BQSR is especially important when    ║
# ║  calling rare variants (AF < 1%) because miscalibrated base qualities    ║
# ║  directly inflate the false positive rate at low-depth positions.        ║
# ║                                                                          ║
# ║  TO ENABLE: Requires known sites VCFs downloaded in Step 0; uncomment.  ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# for SAMPLE in mother father son; do
#     gatk --java-options "-Xmx16g" BaseRecalibrator \
#         -I "${OUTDIR}/${SAMPLE}.markdup.bam" \
#         -R data/ref/Homo_sapiens_assembly38.chr20.fasta \
#         --known-sites data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf.gz \
#         --known-sites data/resources/Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf.gz \
#         --known-sites data/resources/1000G_phase1.snps.high_confidence.hg38.chr20.vcf.gz \
#         -O "${OUTDIR}/${SAMPLE}.recal.table"
#
#     gatk --java-options "-Xmx16g" ApplyBQSR \
#         -I "${OUTDIR}/${SAMPLE}.markdup.bam" \
#         -R data/ref/Homo_sapiens_assembly38.chr20.fasta \
#         --bqsr-recal-file "${OUTDIR}/${SAMPLE}.recal.table" \
#         -O "${OUTDIR}/${SAMPLE}.bqsr.bam"
# done

log ""
log "=== Step 1 complete ==="
log "Outputs in: ${OUTDIR}/"
for SAMPLE in mother father son; do
    if [[ -f "${OUTDIR}/${SAMPLE}.markdup.bam" ]]; then
        log "  ${SAMPLE}.markdup.bam  ✓"
    fi
done
log ""
log "Next step: bash germline/02_haplotypecaller.sh"
