#!/usr/bin/env bash
# =============================================================================
# somatic/04_annotation.sh
# Step 4 of 4 — Somatic Variant Annotation
#
# WHAT THIS STEP DOES:
#   Annotates each PASS somatic variant with functional context from multiple
#   databases, transforming a raw VCF into a clinically interpretable table.
#
# WHY IS ANNOTATION A SEPARATE STEP?
#   Variant calling (Steps 1-3) answers: "Does this position differ from
#   reference in the tumor?" Annotation answers: "So what? What gene is this
#   in? Does it alter protein function? Is it known to drive cancer?"
#
#   A PASS variant from Mutect2 might be:
#     - A frameshift in TP53 (known driver, pathogenic)
#     - A synonymous change in a non-coding region (likely passenger)
#     - A missense in BRCA2 with ClinVar Pathogenic classification
#     - A hotspot in KRAS G12D (COSMIC 45,000+ tumors)
#     - A novel variant in a gene with no known cancer association
#
#   Annotation is what enables triage.
#
# ANNOTATION LAYERS APPLIED:
#
#   1. VEP REST API (Ensembl)
#      The primary annotation engine. For each variant, VEP:
#        - Assigns consequence (SO term): missense_variant, frameshift_variant,
#          stop_gained, splice_donor_variant, synonymous_variant, etc.
#        - Reports IMPACT: HIGH / MODERATE / LOW / MODIFIER
#        - Identifies affected gene + canonical transcript
#        - Generates HGVS notation (HGVSc and HGVSp)
#        - Reports exon/intron position (e.g. "Exon 8 of 11")
#        - Runs SIFT and PolyPhen-2 in silico predictors
#        - Reports conservation scores (phyloP, GERP)
#        - Tags regulatory features (promoter, enhancer, CTCF binding sites)
#
#   2. gnomAD GraphQL API
#      Looks up allele frequency across 125k exomes + 71k genomes.
#      Reports global AF plus 8 ancestry population AFs.
#      A somatic variant with high gnomAD AF (>1%) is likely a germline
#      polymorphism that escaped the normal subtraction — flag for review.
#
#   3. ClinVar eUtils API
#      If the variant has a ClinVar entry, reports:
#        - Clinical significance (Pathogenic/VUS/Benign)
#        - Review status (star rating 1-4)
#        - Associated disease (e.g. "Hereditary breast and ovarian cancer")
#      In somatic calling, ClinVar Pathogenic + PASS from Mutect2 is strong
#      evidence for a clinically relevant variant.
#
#   4. CADD API
#      Combined Annotation Dependent Depletion. A machine learning score
#      integrating 60+ annotation features. CADD phred ≥ 20 = top 1%
#      most deleterious across the genome; ≥ 30 = top 0.1%.
#      Useful for variants not in ClinVar or COSMIC.
#
#   5. CIViC GraphQL API
#      Clinical Interpretation of Variants in Cancer. Community-curated
#      evidence for somatic variants with therapeutic significance.
#      If a variant has CIViC evidence, reports drug associations and
#      evidence level (A = validated, B = clinical trial, C = case study).
#
# INPUTS:  filtered.vcf.gz (PASS somatic variants from Step 3)
# OUTPUTS: annotated.tsv (35-column table), annotated.vcf.gz
#
# RUN:     bash somatic/04_annotation.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

FILTERED_VCF="${ROOT}/results/somatic/filtered/filtered.vcf.gz"
OUTDIR="${ROOT}/results/somatic/annotated"
LOGDIR="${ROOT}/logs/somatic"

mkdir -p "${OUTDIR}" "${LOGDIR}"
log() { printf "[%s] somatic/04 | %s\n" "$(date '+%H:%M:%S')" "$*" | tee -a "${LOGDIR}/04_annotation.log"; }

log "=== Somatic Annotation Pipeline ==="
log "Input: ${FILTERED_VCF}"

NPASS=$(bcftools view -H -f PASS "${FILTERED_VCF}" | wc -l)
log "Variants to annotate: ${NPASS} PASS somatic variants"

# ── Extract PASS variants only ────────────────────────────────────────────────
log ""
log "Extracting PASS variants..."
bcftools view -f PASS -O z -o "${OUTDIR}/pass_only.vcf.gz" "${FILTERED_VCF}"
bcftools index -t "${OUTDIR}/pass_only.vcf.gz"

# ── Run annotation pipeline ───────────────────────────────────────────────────
log ""
log "Running multi-database annotation pipeline..."
log "  VEP REST API   → consequence, HGVS, gene, SIFT/PolyPhen, regulatory"
log "  gnomAD API     → 8-population allele frequencies"
log "  ClinVar eUtils → clinical significance, disease association"
log "  CADD API       → deleteriousness score"
log "  CIViC API      → cancer drug evidence (somatic)"
log ""
log "Note: API calls are rate-limited; annotation takes ~1-2 min per 100 variants"

python "${ROOT}/annotation/annotate_variants.py" \
    --vcf "${OUTDIR}/pass_only.vcf.gz" \
    --type somatic \
    --output "${OUTDIR}/annotated.tsv" \
    --vcf-output "${OUTDIR}/annotated.vcf.gz" \
    --threads 2 \
    2>&1 | tee -a "${LOGDIR}/04_annotation.log"

# ── Annotation summary ─────────────────────────────────────────────────────────
log ""
log "=== Annotation Summary ==="

if [[ -f "${OUTDIR}/annotated.tsv" ]]; then
    NANN=$(tail -n +2 "${OUTDIR}/annotated.tsv" | wc -l)
    log "  Variants annotated: ${NANN}"

    log ""
    log "  Consequence distribution (SO terms):"
    tail -n +2 "${OUTDIR}/annotated.tsv" | \
        awk -F'\t' 'NR>1 && $9!="." {print $9}' | \
        tr "," "\n" | sort | uniq -c | sort -rn | head -10 | \
        awk '{printf "    %-45s %d\n", $2, $1}'

    log ""
    log "  IMPACT distribution:"
    tail -n +2 "${OUTDIR}/annotated.tsv" | \
        awk -F'\t' '$10!="." {print $10}' | \
        sort | uniq -c | sort -rn | \
        awk '{printf "    %-12s %d\n", $2, $1}'

    log ""
    log "  ClinVar clinical significance:"
    tail -n +2 "${OUTDIR}/annotated.tsv" | \
        awk -F'\t' '$25!="." && $25!="" {print $25}' | \
        sort | uniq -c | sort -rn | \
        awk '{printf "    %-35s %d\n", $2, $1}'

    log ""
    log "  Rare variants (gnomAD AF < 0.001):"
    RARE=$(tail -n +2 "${OUTDIR}/annotated.tsv" | \
        awk -F'\t' '$18!="." && $18!="" && $18+0 < 0.001 {count++} END{print count+0}')
    log "    ${RARE} variants with gnomAD AF < 0.1%"
fi

log ""
log "=== Step 4 complete ==="
log "Outputs:"
log "  ${OUTDIR}/annotated.tsv      — full annotation table (35 columns)"
log "  ${OUTDIR}/annotated.vcf.gz   — annotated VCF with INFO fields"
log ""
log "Column reference for annotated.tsv:"
head -1 "${OUTDIR}/annotated.tsv" 2>/dev/null | tr '\t' '\n' | \
    awk '{printf "  [%02d] %s\n", NR, $0}' | \
    tee -a "${LOGDIR}/04_annotation.log"
log ""
log "Next step: jupyter lab somatic/notebooks/somatic_results.ipynb"
log "     -- or: python agent/variant_interpreter.py --tsv ${OUTDIR}/annotated.tsv --type somatic --top 10"
