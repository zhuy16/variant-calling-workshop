#!/usr/bin/env bash
# =============================================================================
# germline/04_hardfilter_annotate.sh
# Step 4 of 4 — Hard Filtering + Annotation
#
# WHAT THIS STEP DOES:
#   1. Separates SNPs and Indels (they have different filter thresholds)
#   2. Applies GATK best-practice hard filters to each variant type
#   3. Merges filtered SNP and Indel VCFs
#   4. Annotates PASS variants with the full multi-database annotation stack
#
# WHY HARD FILTERS INSTEAD OF VQSR?
#   VQSR (Variant Quality Score Recalibration) trains a Gaussian mixture model
#   on variant features (QD, MQ, FS, etc.) using validated truth resources
#   (HapMap, Omni, 1000G). It draws a probabilistic decision boundary between
#   true variants and artifacts.
#
#   VQSR REQUIREMENTS that we don't meet:
#     - Minimum ~10,000 SNPs to train the GMM (our chr20 subset has far fewer)
#     - Truth resource VCFs covering the called regions (we have chr20 subsets
#       but the model needs genome-wide data for robust training)
#     - Typically requires 30+ WGS samples for indel VQSR
#
#   WHEN TO USE HARD FILTERS (this script):
#     - Small panels (targeted sequencing, exome subsets)
#     - Single-sample or small cohort calling (< 30 samples)
#     - Demo data (our case)
#     - Diagnostic labs with small gene panels
#
#   This is NOT a second-best approach — diagnostic labs routinely use hard
#   filters for targeted panels because VQSR is genuinely inappropriate there.
#
# THE HARD FILTER THRESHOLDS AND THEIR MEANING:
#
#   For SNPs:
#   QD < 2.0    QualByDepth: QUAL divided by total depth. A high absolute
#               QUAL on a low-depth site may just reflect deep coverage of
#               an artifact. QD normalizes for depth. QD < 2 = suspicious.
#
#   MQ < 40.0   RMS Mapping Quality: average mapping quality of reads
#               supporting this site. Low MQ means reads map to multiple
#               locations — variants in repetitive regions score low.
#
#   FS > 60.0   FisherStrand: Phred-scaled p-value of strand bias test.
#               Real variants show roughly equal alt support on both strands.
#               FS > 60 = extreme strand bias = likely PCR or mapping artifact.
#
#   SOR > 3.0   StrandOddsRatio: an alternative strand bias metric less
#               sensitive to allele count imbalance than FS. Complements FS.
#
#   MQRankSum < -12.5   Difference in mapping quality between ref and alt reads.
#               If alt reads map with systematically lower MQ than ref reads,
#               the variant may be in a repetitive region affecting only one
#               allele.
#
#   ReadPosRankSum < -8.0   Position of alt alleles within reads. Real variants
#               occur throughout reads; artifacts cluster at read ends where
#               base quality is lowest.
#
#   For Indels (relaxed thresholds — indels are harder to call):
#   QD < 2.0      Same as SNPs
#   FS > 200.0    Much more lenient — short indels commonly show strand bias
#   SOR > 10.0
#   ReadPosRankSum < -20.0
#
# INPUTS:  raw_variants.vcf.gz (from Step 3)
# OUTPUTS: hardfiltered.vcf.gz, annotated.tsv
#
# RUN:     bash germline/04_hardfilter_annotate.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

GENOTYPED="${ROOT}/results/germline/genotyped/raw_variants.vcf.gz"
OUTDIR="${ROOT}/results/germline/filtered"
ANN_DIR="${ROOT}/results/germline/annotated"
LOGDIR="${ROOT}/logs/germline"
REF="${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.fasta"
MEM="-Xmx8g -XX:+UseParallelGC"

mkdir -p "${OUTDIR}" "${ANN_DIR}" "${LOGDIR}"
log() { printf "[%s] germline/04 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/04_hardfilter.log"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)     GENOTYPED="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2";    shift 2 ;;
        *) shift ;;
    esac
done

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SKIPPED: VQSR (Variant Quality Score Recalibration)                    ║
# ║                                                                          ║
# ║  WHAT IT DOES: Trains a Gaussian mixture model on variant annotation    ║
# ║  features using validated truth resources (HapMap 3.3, Omni 2.5,        ║
# ║  1000G high-confidence SNPs). The model learns what "real" variants     ║
# ║  look like in multidimensional feature space (QD, MQ, FS, etc.) and    ║
# ║  assigns each variant a log-odds score (VQSLOD). Variants are then      ║
# ║  filtered by sensitivity tranche (e.g. 99.5% sensitivity tranche        ║
# ║  means we accept the threshold that retains 99.5% of truth variants).   ║
# ║                                                                          ║
# ║  WHY SKIPPED: Our chr20 tutorial callset has ~5,000-10,000 SNPs.         ║
# ║  VQSR needs at minimum ~10,000 SNPs to fit a GMM that generalizes.      ║
# ║  On our data it will either fail with convergence errors or produce a    ║
# ║  poorly fitted model worse than hard filters.                            ║
# ║                                                                          ║
# ║  KEY CONCEPT — Tranches:                                                 ║
# ║    99.9% tranche: very sensitive, more FPs allowed (research use)        ║
# ║    99.0% tranche: balanced (population studies)                          ║
# ║    90.0% tranche: highly specific, lower sensitivity (clinical use)      ║
# ║                                                                          ║
# ║  QUALITY CHECK — Ts/Tv as proxy:                                         ║
# ║    Pre-VQSR: Ts/Tv ≈ 1.9 (some artifacts with random Tv bias)           ║
# ║    Post-VQSR 99.5%: Ts/Tv ≈ 2.05 (artifacts removed, Ti enriched)      ║
# ║                                                                          ║
# ║  TO ENABLE: Requires WGS data with genome-wide truth VCFs.              ║
# ╚══════════════════════════════════════════════════════════════════════════╝
#
# gatk --java-options "${MEM}" VariantRecalibrator \
#     -R "${REF}" \
#     -V "${GENOTYPED}" \
#     --resource:hapmap,known=false,training=true,truth=true,prior=15 hapmap_3.3.hg38.vcf.gz \
#     --resource:omni,known=false,training=true,truth=false,prior=12 1000G_omni2.5.hg38.vcf.gz \
#     --resource:1000G,known=false,training=true,truth=false,prior=10 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#     --resource:dbsnp,known=true,training=false,truth=false,prior=2 dbsnp138.vcf.gz \
#     -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
#     -mode SNP \
#     -O snps.recal \
#     --tranches-file snps.tranches

# ── Step 4A: Separate SNPs and Indels ─────────────────────────────────────────
log "Separating SNPs and Indels for type-specific filtering..."

gatk --java-options "${MEM}" SelectVariants \
    -R "${REF}" \
    -V "${GENOTYPED}" \
    --select-type-to-include SNP \
    -O "${OUTDIR}/raw_snps.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

gatk --java-options "${MEM}" SelectVariants \
    -R "${REF}" \
    -V "${GENOTYPED}" \
    --select-type-to-include INDEL \
    -O "${OUTDIR}/raw_indels.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

log "  SNPs:   $(bcftools view -H "${OUTDIR}/raw_snps.vcf.gz" | wc -l)"
log "  Indels: $(bcftools view -H "${OUTDIR}/raw_indels.vcf.gz" | wc -l)"

# ── Step 4B: Hard Filter — SNPs ───────────────────────────────────────────────
log ""
log "Applying SNP hard filters..."

gatk --java-options "${MEM}" VariantFiltration \
    -R "${REF}" \
    -V "${OUTDIR}/raw_snps.vcf.gz" \
    --filter-expression "QD < 2.0"          --filter-name "QD2" \
    --filter-expression "MQ < 40.0"         --filter-name "MQ40" \
    --filter-expression "FS > 60.0"         --filter-name "FS60" \
    --filter-expression "SOR > 3.0"         --filter-name "SOR3" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${OUTDIR}/filtered_snps.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

SNP_PASS=$(bcftools view -H -f PASS "${OUTDIR}/filtered_snps.vcf.gz" | wc -l)
SNP_FAIL=$(bcftools view -H "${OUTDIR}/filtered_snps.vcf.gz" | wc -l)
log "  SNP PASS: ${SNP_PASS} / ${SNP_FAIL} total"

# ── Step 4C: Hard Filter — Indels ─────────────────────────────────────────────
log ""
log "Applying Indel hard filters (thresholds relaxed vs. SNPs)..."

gatk --java-options "${MEM}" VariantFiltration \
    -R "${REF}" \
    -V "${OUTDIR}/raw_indels.vcf.gz" \
    --filter-expression "QD < 2.0"             --filter-name "QD2" \
    --filter-expression "FS > 200.0"           --filter-name "FS200" \
    --filter-expression "SOR > 10.0"           --filter-name "SOR10" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${OUTDIR}/filtered_indels.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

INDEL_PASS=$(bcftools view -H -f PASS "${OUTDIR}/filtered_indels.vcf.gz" | wc -l)
INDEL_FAIL=$(bcftools view -H "${OUTDIR}/filtered_indels.vcf.gz" | wc -l)
log "  Indel PASS: ${INDEL_PASS} / ${INDEL_FAIL} total"

# ── Step 4D: Merge filtered SNPs and Indels ───────────────────────────────────
log ""
log "Merging filtered SNPs and Indels..."

gatk --java-options "${MEM}" MergeVcfs \
    -I "${OUTDIR}/filtered_snps.vcf.gz" \
    -I "${OUTDIR}/filtered_indels.vcf.gz" \
    -O "${OUTDIR}/hardfiltered.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

TOTAL_PASS=$((SNP_PASS + INDEL_PASS))
log "  Final PASS variants: ${TOTAL_PASS} (${SNP_PASS} SNPs + ${INDEL_PASS} Indels)"

# ── Step 4E: Annotate PASS variants ──────────────────────────────────────────
log ""
log "Running multi-database annotation on PASS germline variants..."

python "${ROOT}/annotation/annotate_variants.py" \
    --vcf "${OUTDIR}/hardfiltered.vcf.gz" \
    --type germline \
    --output "${ANN_DIR}/annotated.tsv" \
    --vcf-output "${ANN_DIR}/annotated.vcf.gz" \
    2>&1 | tee -a "${LOGDIR}/04_hardfilter.log"

# ── Final QC summary ──────────────────────────────────────────────────────────
log ""
log "=== Final Germline Callset Summary ==="
log "  PASS SNPs:    ${SNP_PASS}"
log "  PASS Indels:  ${INDEL_PASS}"
log "  Total PASS:   ${TOTAL_PASS}"

if [[ -f "${ANN_DIR}/annotated.tsv" ]]; then
    log ""
    log "  Annotation summary:"
    log "  ClinVar significance distribution:"
    tail -n +2 "${ANN_DIR}/annotated.tsv" | \
        awk -F'\t' '$25!="." && $25!="" {print $25}' | \
        sort | uniq -c | sort -rn | \
        awk '{printf "    %-35s %d\n", $2, $1}'

    RARE=$(tail -n +2 "${ANN_DIR}/annotated.tsv" | \
        awk -F'\t' '$18!="." && $18!="" && $18+0 < 0.001 {c++} END{print c+0}')
    log ""
    log "  Rare variants (gnomAD AF < 0.1%): ${RARE}"
    log "  These are candidates for rare disease analysis"
fi

log ""
log "=== Step 4 complete ==="
log "Outputs:"
log "  ${OUTDIR}/hardfiltered.vcf.gz     — PASS germline variants"
log "  ${ANN_DIR}/annotated.tsv          — fully annotated table"
log ""
log "Next steps:"
log "  jupyter lab germline/notebooks/germline_results.ipynb"
log "  python agent/variant_interpreter.py --tsv ${ANN_DIR}/annotated.tsv --type germline --top 10"
