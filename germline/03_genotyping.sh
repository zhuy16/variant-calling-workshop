#!/usr/bin/env bash
# =============================================================================
# germline/03_genotyping.sh
# Step 3 of 4 — Joint Genotyping (CombineGVCFs + GenotypeGVCFs)
#
# WHAT THIS STEP DOES:
#   1. Merges the three per-sample GVCFs into a single cohort GVCF
#   2. Performs joint genotyping to produce a multi-sample VCF
#
# WHY JOINT GENOTYPING MATTERS:
#   Consider a rare variant at chr20:500000 present only in the son.
#
#   SINGLE-SAMPLE calling (son only):
#     Son has 3 reads supporting the alt allele out of 20 total (15% VAF).
#     HaplotypeCaller alone is uncertain — could be a real rare het or an error.
#     Likely filtered as low-quality or not called.
#
#   JOINT GENOTYPING (mother + father + son):
#     Mother: 0 alt reads out of 25 → confident reference (GQ=50)
#     Father: 0 alt reads out of 28 → confident reference (GQ=55)
#     Son:    3 alt reads out of 20 → uncertain (GQ=18)
#
#     Now the joint model knows:
#       - Both parents are confidently reference
#       - The alt allele has allele count 1 in this cohort (very rare)
#       - The son's genotype, though low GQ, is the only explanation for
#         the observed alt reads given the parental genotypes
#       → The variant is called with higher confidence as a de novo candidate
#
#   This "borrowing strength" across samples is the core advantage of the
#   GVCF/joint genotyping workflow over per-sample calling.
#
# COMBINEAGVCFS vs. GENOMICSDBIMPORT:
#   For small cohorts (< 100 samples), CombineGVCFs is simpler.
#   For large cohorts (hundreds to thousands), GenomicsDBImport is required
#   as it uses a columnar database backend that scales much better.
#   In production genomics (UK Biobank, gnomAD), GenomicsDBImport is standard.
#
# GENOTYPEGVCFS — what it outputs per variant:
#   CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  MOTHER  FATHER  SON
#   chr20  500000  .  A  T  42  .  AC=1;AF=0.167;AN=6  GT:AD:GQ:PL
#         0/0:25,0:50:0,75,900   0/0:28,0:55:0,84,990   0/1:17,3:18:54,0,450
#
#   INFO fields added by GenotypeGVCFs:
#     AC   — Allele count in genotypes (1 = only one copy across all samples)
#     AF   — Allele frequency (1/6 = 0.167 in a trio)
#     AN   — Total alleles genotyped (6 = 2 per sample × 3 samples)
#     DP   — Combined depth across all samples
#     QD   — QualByDepth (QUAL / DP) — key filter metric
#     MQ   — RMS mapping quality of reads
#     FS   — FisherStrand (strand bias test)
#     SOR  — StrandOddsRatio (alternative strand bias test)
#     MQRankSum — Mapping quality difference between ref and alt reads
#     ReadPosRankSum — Position of alt alleles within reads vs. ref reads
#
# INPUTS:  {mother,father,son}.g.vcf.gz (per-sample GVCFs from Step 2)
# OUTPUTS: raw_variants.vcf.gz (joint-genotyped, pre-filter multi-sample VCF)
#
# RUN:     bash germline/03_genotyping.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

GVCF_DIR="${ROOT}/results/germline/gvcfs"
OUTDIR="${ROOT}/results/germline/genotyped"
LOGDIR="${ROOT}/logs/germline"
REF="${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.fasta"
CHROM="20"   # b37 reference naming (SN:20)
DBSNP="${ROOT}/data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf"
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] germline/03 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/03_genotyping.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --gvcf-args) GVCF_ARGS="$2"; shift 2 ;;
        --outdir)    OUTDIR="$2";    shift 2 ;;
        *) shift ;;
    esac
done

# ── Build GVCF input args ─────────────────────────────────────────────────────
if [[ -z "${GVCF_ARGS:-}" ]]; then
    GVCF_ARGS=""
    for SAMPLE in mother father son; do
        GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
        if [[ ! -f "${GVCF}" ]]; then
            log "ERROR: ${GVCF} not found. Run Step 2 first."
            exit 1
        fi
        GVCF_ARGS="${GVCF_ARGS} -V ${GVCF}"
    done
fi

# ── Step 3A: CombineGVCFs ────────────────────────────────────────────────────
log "Merging per-sample GVCFs into cohort GVCF..."
log "  Samples: mother (NA12878), father (NA12891), son (NA12892)"
log "  This preserves all reference confidence blocks needed for joint calling"

gatk --java-options "${MEM}" CombineGVCFs \
    -R "${REF}" \
    ${GVCF_ARGS} \
    -L "${CHROM}" \
    -O "${OUTDIR}/cohort.g.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/03_genotyping.log"

log "Cohort GVCF created."
COHORT_SIZE=$(du -sh "${OUTDIR}/cohort.g.vcf.gz" | cut -f1)
log "  Size: ${COHORT_SIZE}"

# ── Step 3B: GenotypeGVCFs ────────────────────────────────────────────────────
log ""
log "Running joint genotyping (GenotypeGVCFs)..."
log "  This is where shared statistical power across samples is applied"
log "  Output will be a multi-sample VCF: one column per trio member"

gatk --java-options "${MEM}" GenotypeGVCFs \
    -R "${REF}" \
    -V "${OUTDIR}/cohort.g.vcf.gz" \
    --dbsnp "${DBSNP}" \
    -L "${CHROM}" \
    -O "${OUTDIR}/raw_variants.vcf.gz" \
    --standard-min-confidence-threshold-for-calling 30 \
    2>&1 | tee -a "${LOGDIR}/03_genotyping.log"

# ── QC: Callset statistics ────────────────────────────────────────────────────
log ""
log "=== Raw Callset Statistics ==="

TOTAL=$(bcftools view -H "${OUTDIR}/raw_variants.vcf.gz" | wc -l)
SNPS=$(bcftools view -H -v snps "${OUTDIR}/raw_variants.vcf.gz" | wc -l)
INDELS=$(bcftools view -H -v indels "${OUTDIR}/raw_variants.vcf.gz" | wc -l)

log "  Total variants:  ${TOTAL}"
log "  SNPs:            ${SNPS}"
log "  Indels:          ${INDELS}"
log "  SNP/Indel ratio: $(awk "BEGIN{printf \"%.1f\", ${SNPS}/${INDELS:-1}}")"

# Compute Ti/Tv ratio — a key callset quality metric
log ""
log "  Computing Ts/Tv (transition/transversion) ratio..."
log "  Expected for WGS: ~2.0-2.1 | For WES: ~2.6-3.5"
log "  Low Ts/Tv (<1.8) suggests excess false positives"
bcftools stats "${OUTDIR}/raw_variants.vcf.gz" 2>/dev/null | \
    awk '/^TSTV/{printf "  Ts/Tv = %.2f  (ts=%s, tv=%s)\n", $5, $3, $4}' || true

log ""
log "  Allele count distribution (AC):"
log "  AC=1: rare/private variants (present in 1/6 alleles = only one het sample)"
log "  AC=2: present in 2 alleles (one hom or two het samples)"
log "  AC=6: all samples homozygous alt"
bcftools query -f '%INFO/AC\n' "${OUTDIR}/raw_variants.vcf.gz" | \
    sort | uniq -c | sort -rn | head -5 | \
    awk '{printf "  AC=%-4s count=%d\n", $2, $1}'

log ""
log "  INTERPRETING RAW CALLSET:"
log "  The raw callset contains real variants + artifacts."
log "  Artifacts include: PCR errors, alignment errors near repeats, strand bias."
log "  Step 4 applies hard filters to remove these — see 04_hardfilter_annotate.sh"

log ""
log "=== Step 3 complete ==="
log "Outputs:"
log "  ${OUTDIR}/cohort.g.vcf.gz    — combined per-sample GVCFs"
log "  ${OUTDIR}/raw_variants.vcf.gz — joint-genotyped variants (pre-filter)"
log ""
log "Next step: bash germline/04_hardfilter_annotate.sh"
