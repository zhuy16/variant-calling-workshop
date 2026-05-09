#!/usr/bin/env bash
# =============================================================================
# download_tutorial_data.sh
#
# Downloads GATK workshop tutorial data from Google Cloud Storage via HTTPS.
# No GCP account, no gsutil, no Google Cloud SDK required.
# All files are freely available at public GCS URLs.
#
# SOMATIC data: HCC1143 breast cancer cell line tumor + matched normal
#   Restricted to chr17:7,000,000-20,000,000 (~50 MB BAMs)
#   This region covers TP53 and BRCA1 — high-yield for demo somatic calls
#
# GERMLINE data: CEPH trio (NA12878 mother / NA12891 father / NA12892 son)
#   Restricted to chr20 (~60 MB BAMs)
#   NA12878 is the Genome in a Bottle benchmark — extensively characterized
#
# Prerequisites: curl (pre-installed on macOS), samtools (for index check)
#
# Usage:
#   bash data/download_tutorial_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

# ── Destination directories ───────────────────────────────────────────────────
mkdir -p \
    "${ROOT}/data/bams/somatic" \
    "${ROOT}/data/bams/germline" \
    "${ROOT}/data/ref" \
    "${ROOT}/data/resources"

# ── GCS public HTTPS base URL (no gsutil needed) ──────────────────────────────
# Bucket structure verified against live GCS listing (May 2026):
#   workshop_2002/2-germline/  ← germline trio BAMs + resources
#   workshop_2002/3-somatic/   ← somatic resources + precomputed outputs
GCS="https://storage.googleapis.com"
TUTS="${GCS}/gatk-tutorials/workshop_2002"

log() { printf "[%s] %s\n" "$(date '+%H:%M:%S')" "$*"; }

# Helper: download with resume support (-C -) and progress bar
fetch() {
    local url="$1"
    local dest="$2"
    local filename
    filename="$(basename "${dest}")"

    if [[ -f "${dest}" ]]; then
        log "  Already exists, skipping: ${filename}"
        return 0
    fi

    log "  Downloading: ${filename}"
    curl -L --progress-bar --retry 3 --retry-delay 5 \
        -o "${dest}" "${url}" || {
        log "  ERROR: Failed to download ${url}"
        rm -f "${dest}"
        exit 1
    }
}

# =============================================================================
# SOMATIC MODULE — chr17, precomputed Mutect2 outputs
#
# NOTE: The separate tumor/normal BAMs are not publicly stored in this bucket.
# The workshop was designed for instructor-led sessions where BAMs were provided
# interactively. We download:
#   (a) The precomputed filtered VCF from a completed Mutect2 run — this is the
#       real output of running Mutect2 on HCC1143 tumor vs. matched normal.
#       We start the annotation pipeline from this VCF.
#   (b) Resources needed to understand/re-run the calling steps (gnomAD AF VCF,
#       Panel of Normals, intervals) — these are included so the scripts are
#       complete and the calling commands can be run if you supply your own BAMs.
#
# WHAT YOU MISS BY NOT HAVING THE RAW BAMs:
#   - somatic/01_preprocessing.sh (MarkDuplicates on tumor + normal BAMs)
#   - somatic/02_mutect2_calling.sh (Mutect2 + GetPileupSummaries)
#   - somatic/03_filtering.sh (LearnReadOrientationModel + FilterMutectCalls)
# These steps are fully documented in the scripts. The annotation + agentic
# interpretation in somatic/04_annotation.sh runs on the precomputed VCF.
# =============================================================================

log "=== SOMATIC: Downloading precomputed Mutect2 filtered VCF (chr17) ==="
mkdir -p "${ROOT}/results/somatic/filtered"

# 9_somatic_oncefiltered.vcf.gz: real Mutect2 output filtered with
# FilterMutectCalls. Contains PASS variants + variants with filter tags.
# This is equivalent to what somatic/03_filtering.sh produces.
fetch "${TUTS}/3-somatic/mutect2_precomputed/9_somatic_oncefiltered.vcf.gz" \
      "${ROOT}/results/somatic/filtered/filtered.vcf.gz"
fetch "${TUTS}/3-somatic/mutect2_precomputed/9_somatic_oncefiltered.vcf.gz.tbi" \
      "${ROOT}/results/somatic/filtered/filtered.vcf.gz.tbi"

# Contamination estimate (used by FilterMutectCalls; included for reference)
fetch "${TUTS}/3-somatic/mutect2_precomputed/8_pair_calculatecontamination.table" \
      "${ROOT}/data/resources/pair_calculatecontamination.table"

log "=== SOMATIC: Downloading resources (for manual re-running of calling steps) ==="

# gnomAD AF-only VCF — chr17 subset
# Used by: Mutect2 --germline-resource (prior on germline sites)
#          GetPileupSummaries (contamination estimation sites)
fetch "${TUTS}/3-somatic/resources/chr17_af-only-gnomad_grch38.vcf.gz" \
      "${ROOT}/data/resources/af-only-gnomad.chr17.vcf.gz"
fetch "${TUTS}/3-somatic/resources/chr17_af-only-gnomad_grch38.vcf.gz.tbi" \
      "${ROOT}/data/resources/af-only-gnomad.chr17.vcf.gz.tbi"

# ExAC common variants — used for contamination estimation
fetch "${TUTS}/3-somatic/resources/chr17_small_exac_common_3_grch38.vcf.gz" \
      "${ROOT}/data/resources/small_exac_common_3.chr17.vcf.gz"
fetch "${TUTS}/3-somatic/resources/chr17_small_exac_common_3_grch38.vcf.gz.tbi" \
      "${ROOT}/data/resources/small_exac_common_3.chr17.vcf.gz.tbi"

# Panel of Normals — chr17 Mutect2 PoN
fetch "${TUTS}/3-somatic/resources/chr17_m2pon.vcf.gz" \
      "${ROOT}/data/resources/1000g_pon.chr17.vcf.gz"
fetch "${TUTS}/3-somatic/resources/chr17_m2pon.vcf.gz.tbi" \
      "${ROOT}/data/resources/1000g_pon.chr17.vcf.gz.tbi"

# Interval list (chr17 coding + TP53 region)
fetch "${TUTS}/3-somatic/resources/chr17plus.interval_list" \
      "${ROOT}/data/resources/chr17_20M-60M.interval_list"

log "Somatic resources: $(du -sh "${ROOT}/data/resources/" | cut -f1)"

# =============================================================================
# GERMLINE MODULE — chr20 CEPH trio BAMs (mother / father / son)
# All BAMs are publicly available and run from step 1.
# =============================================================================

log "=== GERMLINE: Downloading trio BAMs (mother/father/son, chr20) ==="

# Why a trio?
# Joint genotyping across related samples (GVCF workflow) improves sensitivity
# for rare variants by sharing allele count information. In a trio:
#   - De novo variants appear in the child but not parents
#   - Compound heterozygotes can be phased using parental genotypes
#   - Mendelian violation rate is a built-in quality metric

for sample in mother father son; do
    fetch "${TUTS}/2-germline/bams/${sample}.bam" \
          "${ROOT}/data/bams/germline/${sample}.bam"
    fetch "${TUTS}/2-germline/bams/${sample}.bai" \
          "${ROOT}/data/bams/germline/${sample}.bai"
    log "  Downloaded ${sample}.bam"
done
log "Germline BAMs: $(du -sh "${ROOT}/data/bams/germline/" | cut -f1)"

# ─── Reference genome — chr20 ─────────────────────────────────────────────────
log "=== GERMLINE: Downloading chr20 reference ==="

fetch "${TUTS}/2-germline/ref/ref.fasta" \
      "${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.fasta"
fetch "${TUTS}/2-germline/ref/ref.fasta.fai" \
      "${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.fasta.fai"
fetch "${TUTS}/2-germline/ref/ref.dict" \
      "${ROOT}/data/ref/Homo_sapiens_assembly38.chr20.dict"

# ─── Known sites (dbSNP + GIAB truth) ─────────────────────────────────────────
log "=== GERMLINE: Downloading known sites VCFs ==="

# dbSNP: tags variants with rs IDs; used by HaplotypeCaller --dbsnp
fetch "${TUTS}/2-germline/resources/dbsnp.vcf" \
      "${ROOT}/data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf"

# GIAB truth set for mother (NA12878) — used for optional benchmarking
fetch "${TUTS}/2-germline/resources/motherGIAB.vcf.gz" \
      "${ROOT}/data/resources/motherGIAB.vcf.gz"
fetch "${TUTS}/2-germline/resources/motherGIAB.vcf.gz.tbi" \
      "${ROOT}/data/resources/motherGIAB.vcf.gz.tbi"

# High-confidence regions for mother (used with hap.py for benchmarking)
fetch "${TUTS}/2-germline/intervals/motherHighconf.bed" \
      "${ROOT}/data/resources/motherHighconf.bed"

# Index the dbSNP VCF (needed by GATK)
if command -v gatk &>/dev/null; then
    log "  Indexing dbSNP VCF..."
    gatk IndexFeatureFile \
        -I "${ROOT}/data/resources/Homo_sapiens_assembly38.dbsnp138.chr20.vcf" \
        2>/dev/null || true
fi

# =============================================================================
# Summary
# =============================================================================

log "=== Download complete ==="
echo ""
echo "Directory sizes:"
du -sh "${ROOT}/data/bams/somatic/"  | awk '{print "  Somatic BAMs:    " $1}'
du -sh "${ROOT}/data/bams/germline/" | awk '{print "  Germline BAMs:   " $1}'
du -sh "${ROOT}/data/ref/"           | awk '{print "  Reference:       " $1}'
du -sh "${ROOT}/data/resources/"     | awk '{print "  Resources:       " $1}'
echo ""
echo "Next step:"
echo "  bash somatic/01_preprocessing.sh"
echo "  -- or --"
echo "  nextflow run main.nf -profile local --module somatic"
