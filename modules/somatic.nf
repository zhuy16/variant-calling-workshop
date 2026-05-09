#!/usr/bin/env nextflow
// =============================================================================
// modules/somatic.nf  —  Nextflow DSL2 somatic variant calling module
//
// Wraps the 4-step GATK Mutect2 somatic pipeline as named Nextflow processes.
// Each process corresponds to one bash script in somatic/, making it easy to
// run individual steps for learning while the Nextflow workflow chains them.
// =============================================================================

nextflow.enable.dsl = 2

// ── Load config ───────────────────────────────────────────────────────────────
def cfg = new org.yaml.snakeyaml.Yaml().load(
    new File(params.config).text
)

// ── Process 1: Mark Duplicates ────────────────────────────────────────────────
// Removes PCR and optical duplicate reads before variant calling.
// Without this step, PCR duplicates would inflate apparent depth and distort
// tumor variant allele frequencies (VAF), making clonal variants appear higher
// than they truly are and subclonal variants appear amplified.
process MARK_DUPLICATES_SOMATIC {
    label 'gatk'
    tag { sample }

    publishDir "${params.outdir}/somatic/preprocessing", mode: 'copy', pattern: '*.metrics'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.markdup.bam"), path("${sample}.markdup.bam.bai"), emit: bam
    path "${sample}.markdup.metrics", emit: metrics

    script:
    """
    bash ${projectDir}/somatic/01_preprocessing.sh \\
        --bam ${bam} \\
        --sample ${sample} \\
        --outdir . \\
        --threads ${task.cpus}
    """
}

// ── Process 2: Mutect2 Tumor-Normal Calling ───────────────────────────────────
// Runs Mutect2 in tumor-normal mode, simultaneously calling variants and
// collecting read orientation data (f1r2) needed for the artifact model.
// Also runs GetPileupSummaries + CalculateContamination to quantify
// cross-sample contamination — a critical somatic QC metric.
process MUTECT2_CALL {
    label 'gatk'

    publishDir "${params.outdir}/somatic/mutect2", mode: 'copy'

    input:
    tuple val(tumor_sample), path(tumor_bam), path(tumor_bai)
    tuple val(normal_sample), path(normal_bam), path(normal_bai)

    output:
    path "unfiltered.vcf.gz",           emit: vcf
    path "unfiltered.vcf.gz.tbi",       emit: tbi
    path "f1r2.tar.gz",                 emit: f1r2
    path "tumor_pileups.table",         emit: tumor_pileups
    path "normal_pileups.table",        emit: normal_pileups
    path "contamination.table",         emit: contamination
    path "segments.table",              emit: segments

    script:
    """
    bash ${projectDir}/somatic/02_mutect2_calling.sh \\
        --tumor-bam ${tumor_bam} \\
        --normal-bam ${normal_bam} \\
        --tumor-sample ${tumor_sample} \\
        --normal-sample ${normal_sample} \\
        --outdir . \\
        --threads ${task.cpus}
    """
}

// ── Process 3: Filter Mutect Calls ───────────────────────────────────────────
// Applies the full cascade of Mutect2 filters using:
//   - Contamination table (from CalculateContamination)
//   - Tumor segmentation (local copy number context)
//   - Orientation bias priors (from LearnReadOrientationModel — detects
//     FFPE-associated C→T oxidation artifacts in a strand-specific manner)
process FILTER_MUTECT_CALLS {
    label 'gatk'

    publishDir "${params.outdir}/somatic/filtered", mode: 'copy'

    input:
    path vcf
    path tbi
    path f1r2
    path contamination
    path segments

    output:
    path "filtered.vcf.gz",     emit: vcf
    path "filtered.vcf.gz.tbi", emit: tbi
    path "filtering_stats.tsv", emit: stats

    script:
    """
    bash ${projectDir}/somatic/03_filtering.sh \\
        --vcf ${vcf} \\
        --f1r2 ${f1r2} \\
        --contamination ${contamination} \\
        --segments ${segments} \\
        --outdir .
    """
}

// ── Process 4: Annotate Somatic Variants ─────────────────────────────────────
// Chains the full annotation stack:
//   VEP REST → gnomAD GraphQL → ClinVar eUtils → CADD API → CIViC GraphQL
// Outputs enriched TSV (one variant per row, ~35 annotation columns)
// plus a VCF with INFO fields populated for downstream tools.
process ANNOTATE_SOMATIC {
    label 'python'

    publishDir "${params.outdir}/somatic/annotated", mode: 'copy'

    input:
    path vcf

    output:
    path "annotated.tsv", emit: tsv
    path "annotated.vcf.gz", emit: vcf

    script:
    """
    python ${projectDir}/annotation/annotate_variants.py \\
        --vcf ${vcf} \\
        --type somatic \\
        --output annotated.tsv \\
        --vcf-output annotated.vcf.gz
    """
}

// ── Somatic workflow ──────────────────────────────────────────────────────────
workflow SOMATIC_WORKFLOW {

    def tumor_bam  = file(cfg.paths.somatic.bam_tumor)
    def normal_bam = file(cfg.paths.somatic.bam_normal)
    def tumor_name  = cfg.samples.somatic.tumor_name
    def normal_name = cfg.samples.somatic.normal_name

    tumor_ch  = Channel.of(tuple(tumor_name,  tumor_bam,  file(tumor_bam + ".bai")))
    normal_ch = Channel.of(tuple(normal_name, normal_bam, file(normal_bam + ".bai")))

    MARK_DUPLICATES_SOMATIC(tumor_ch.mix(normal_ch))

    tumor_markdup  = MARK_DUPLICATES_SOMATIC.out.bam.filter { it[0] == tumor_name }
    normal_markdup = MARK_DUPLICATES_SOMATIC.out.bam.filter { it[0] == normal_name }

    MUTECT2_CALL(tumor_markdup, normal_markdup)

    FILTER_MUTECT_CALLS(
        MUTECT2_CALL.out.vcf,
        MUTECT2_CALL.out.tbi,
        MUTECT2_CALL.out.f1r2,
        MUTECT2_CALL.out.contamination,
        MUTECT2_CALL.out.segments
    )

    ANNOTATE_SOMATIC(FILTER_MUTECT_CALLS.out.vcf)
}
