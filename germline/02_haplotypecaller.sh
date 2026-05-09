#!/usr/bin/env bash
# =============================================================================
# germline/02_haplotypecaller.sh
# Step 2 of 4 — HaplotypeCaller in GVCF mode (per sample)
#
# WHAT THIS STEP DOES:
#   Runs HaplotypeCaller on each trio member independently, producing a GVCF
#   (genomic VCF) per sample. These GVCFs are inputs to joint genotyping in
#   Step 3.
#
# WHY GVCF MODE (-ERC GVCF)?
#   Without GVCF mode, HaplotypeCaller produces a regular VCF containing ONLY
#   variant sites. This creates a critical problem for multi-sample analysis:
#
#   Regular VCF (per-sample calling):
#     Sample A: chr20  100  A  T  (variant)
#     Sample B: [chr20 100 not reported — was it reference or not called?]
#     AMBIGUOUS: You cannot distinguish "reference" from "no data"
#
#   GVCF mode:
#     Sample A: chr20  100  A  T  (variant)
#     Sample B: chr20  100  A  <NON_REF>  GQ=45  DP=32  (confident reference)
#     UNAMBIGUOUS: Sample B is definitely reference with high confidence
#
#   This distinction is essential for:
#     - Accurately computing allele counts across a cohort
#     - Detecting de novo variants (requires proof the parent is reference)
#     - Joint genotyping across samples added later (no reprocessing)
#
# HOW HAPLOTYPECALLER WORKS:
#   1. ACTIVE REGION DETECTION:
#      Scans the BAM for windows where the read data deviates significantly
#      from the reference (elevated mismatch rate). These become "active regions"
#      where variant calling is attempted.
#
#   2. LOCAL DE NOVO ASSEMBLY:
#      Within each active region (~200bp window), builds a De Bruijn graph of
#      all observed k-mer sequences. Traversing this graph generates candidate
#      haplotypes — possible sequences that could explain the reads.
#
#   3. PAIRHMM LIKELIHOOD CALCULATION:
#      For each read, calculates the likelihood it was generated from each
#      candidate haplotype using a Pair Hidden Markov Model that accounts for
#      base quality, mapping quality, and indel probabilities.
#
#   4. GENOTYPE LIKELIHOOD ESTIMATION:
#      Marginalizes over haplotype pairs to compute P(genotype | reads):
#        GL(0/0) — probability of homozygous reference
#        GL(0/1) — probability of heterozygous
#        GL(1/1) — probability of homozygous alternate
#      These genotype likelihoods (GL/PL fields) are stored in the GVCF.
#
#   5. GVCF EMISSION:
#      Emits all positions — variant sites with allele-specific likelihoods,
#      and non-variant blocks as reference confidence bands (RR records).
#
# KEY GVCF OUTPUT FIELDS:
#   GT    Genotype (0/1 het, 1/1 hom-alt, ./. no-call)
#   GQ    Genotype quality (Phred-scaled confidence in GT)
#   DP    Total read depth at this site
#   AD    Allele-specific depths [ref_depth, alt_depth]
#   PL    Phred-scaled genotype likelihoods [PL(0/0), PL(0/1), PL(1/1)]
#
# INTERPRETING GQ:
#   GQ ≥ 20  → 99% confident in the genotype (most variant callers require this)
#   GQ ≥ 30  → 99.9% confident (used in clinical contexts)
#   GQ < 10  → low confidence; treat as uncertain; do not report
#
# INPUTS:  {mother,father,son}.markdup.bam (from Step 1)
# OUTPUTS: {mother,father,son}.g.vcf.gz (per-sample GVCFs)
#
# RUN:     bash germline/02_haplotypecaller.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

PREPDIR="${ROOT}/results/germline/preprocessing"
OUTDIR="${ROOT}/results/germline/gvcfs"
LOGDIR="${ROOT}/logs/germline"
REF="${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.fasta"
CHROM="20"   # reference uses b37 naming (SN:20, not chr20)
DBSNP="${ROOT}/data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf"
THREADS=4
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] germline/02 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/02_haplotypecaller.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)     INPUT_BAM="$2"; shift 2 ;;
        --sample)  SAMPLE="$2";    shift 2 ;;
        --outdir)  OUTDIR="$2";    shift 2 ;;
        --threads) THREADS="$2";   shift 2 ;;
        *) shift ;;
    esac
done

run_haplotypecaller() {
    local SAMPLE="$1"
    local BAM="$2"

    log "Running HaplotypeCaller on ${SAMPLE}..."
    log "  BAM: ${BAM}"
    log "  Output: ${OUTDIR}/${SAMPLE}.g.vcf.gz"
    log "  Mode: -ERC GVCF (emits reference confidence records)"

    gatk --java-options "${MEM}" HaplotypeCaller \
        -R "${REF}" \
        -I "${BAM}" \
        -O "${OUTDIR}/${SAMPLE}.g.vcf.gz" \
        -ERC GVCF \
        -L "${CHROM}" \
        --dbsnp "${DBSNP}" \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 30 \
        --native-pair-hmm-threads "${THREADS}" \
        2>&1 | tee -a "${LOGDIR}/02_haplotypecaller.log"

    local NVAR=$(bcftools view -H -e 'ALT=="<NON_REF>"' "${OUTDIR}/${SAMPLE}.g.vcf.gz" | wc -l)
    log "  ${SAMPLE}: ${NVAR} variant sites called (plus reference blocks)"
}

# ── If called directly, run all three samples sequentially ───────────────────
# (In Nextflow, these run in parallel — one process per sample)
if [[ -z "${INPUT_BAM:-}" ]]; then
    log "Running HaplotypeCaller on all trio members..."
    log "NOTE: In Nextflow, these run in parallel. Running sequentially here."
    log ""

    for SAMPLE in mother father son; do
        BAM="${PREPDIR}/${SAMPLE}.markdup.bam"
        if [[ ! -f "${BAM}" ]]; then
            log "ERROR: ${BAM} not found. Run Step 1 first."
            exit 1
        fi
        run_haplotypecaller "${SAMPLE}" "${BAM}"
        log ""
    done
else
    run_haplotypecaller "${SAMPLE:-unknown}" "${INPUT_BAM}"
fi

# ── QC: Inspect one GVCF to understand format ─────────────────────────────────
log "=== GVCF Format Inspection (mother, first 5 variant records) ==="
MOTHER_GVCF="${OUTDIR}/mother.g.vcf.gz"
if [[ -f "${MOTHER_GVCF}" ]]; then
    log "First 5 variant sites in mother.g.vcf.gz:"
    bcftools view -H -e 'ALT=="<NON_REF>"' "${MOTHER_GVCF}" | head -5 | \
        awk '{printf "  chr20:%s %s>%s GT=%s GQ=%s DP=%s\n",
              $2, $4, $5,
              gensub(/.*GT:/, "", 1, $9) == "" ? "?" : substr($10, 1, 3),
              gensub(/.*:/, "", 1, gensub(/:([^:]+):.*/, "\\1", 1, $10)),
              gensub(/.*:.*:([^:]+):.*/, "\\1", 1, $10)}' \
        2>/dev/null || \
    bcftools view -H -e 'ALT=="<NON_REF>"' "${MOTHER_GVCF}" | head -5 | \
        awk '{print "  " $1 ":" $2 " " $4 ">" $5 " FILTER=" $7}'

    log ""
    log "  Field guide:"
    log "    ALT=<NON_REF>         → reference confidence block (non-variant)"
    log "    ALT=A,<NON_REF>       → variant site; <NON_REF> is a symbolic allele"
    log "    GT=0/1                → heterozygous variant"
    log "    GQ=45                 → 99.997% confident this is a het, not ref or hom"
    log "    AD=30,15              → 30 ref reads, 15 alt reads (50% VAF = expected het)"
    log "    PL=450,0,900          → P(ref)=10^-45, P(het)=10^0=1, P(hom)=10^-90"
fi

log ""
log "=== Step 2 complete ==="
log "Outputs:"
for SAMPLE in mother father son; do
    if [[ -f "${OUTDIR}/${SAMPLE}.g.vcf.gz" ]]; then
        SIZE=$(du -sh "${OUTDIR}/${SAMPLE}.g.vcf.gz" | cut -f1)
        log "  ${SAMPLE}.g.vcf.gz  (${SIZE})"
    fi
done
log ""
log "Next step: bash germline/03_genotyping.sh"
