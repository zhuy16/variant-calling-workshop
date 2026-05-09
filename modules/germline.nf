#!/usr/bin/env nextflow
// =============================================================================
// modules/germline.nf  —  Nextflow DSL2 germline variant calling module
//
// Implements the GATK GVCF workflow (HaplotypeCaller → CombineGVCFs →
// GenotypeGVCFs → hard filters → annotation) for a mother/father/son trio.
//
// A key advantage of Nextflow here: the 3 HaplotypeCaller jobs (one per sample)
// run in PARALLEL automatically via process { } channel semantics.
// In bash you'd need to manage background jobs manually.
// =============================================================================

nextflow.enable.dsl = 2

def cfg = new org.yaml.snakeyaml.Yaml().load(
    new File(params.config).text
)

// ── Process 1: Mark Duplicates (per sample) ───────────────────────────────────
// Identical in concept to the somatic preprocessing step. Run on all 3 trio
// members independently — Nextflow parallelises these automatically.
process MARK_DUPLICATES_GERMLINE {
    label 'gatk'
    tag { sample }

    publishDir "${params.outdir}/germline/preprocessing/${sample}", mode: 'copy', pattern: '*.metrics'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.markdup.bam"), path("${sample}.markdup.bam.bai"), emit: bam
    path "${sample}.markdup.metrics", emit: metrics

    script:
    """
    bash ${projectDir}/germline/01_preprocessing.sh \\
        --bam ${bam} \\
        --sample ${sample} \\
        --outdir . \\
        --threads ${task.cpus}
    """
}

// ── Process 2: HaplotypeCaller GVCF mode (per sample) ────────────────────────
// Generates a GVCF (genomic VCF) for each sample individually.
// GVCF mode (-ERC GVCF) emits reference confidence records for all positions,
// not just variant sites. This is what enables scalable joint genotyping:
//   - New samples can be added without reprocessing existing ones
//   - Statistical power is shared across samples at joint genotyping step
//   - Variant calling model uses local de novo assembly (active regions)
process HAPLOTYPECALLER_GVCF {
    label 'gatk'
    tag { sample }

    publishDir "${params.outdir}/germline/gvcfs", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi")

    script:
    """
    bash ${projectDir}/germline/02_haplotypecaller.sh \\
        --bam ${bam} \\
        --sample ${sample} \\
        --outdir . \\
        --threads ${task.cpus}
    """
}

// ── Process 3: Joint Genotyping (CombineGVCFs + GenotypeGVCFs) ───────────────
// Combines per-sample GVCFs and performs joint genotyping.
//
// Why joint genotyping matters:
//   If a variant is present in only one sample at low coverage, calling it
//   from that sample alone is uncertain. Joint genotyping sees that the other
//   two samples have reference alleles at adequate depth, which increases
//   confidence the variant is real (not a sequencing artifact) or confirms
//   it as a de novo if parents are truly homozygous reference.
process JOINT_GENOTYPING {
    label 'gatk'

    publishDir "${params.outdir}/germline/genotyped", mode: 'copy'

    input:
    path gvcfs
    path tbis

    output:
    path "raw_variants.vcf.gz",     emit: vcf
    path "raw_variants.vcf.gz.tbi", emit: tbi

    script:
    def gvcf_args = gvcfs.collect { "-V ${it}" }.join(' ')
    """
    bash ${projectDir}/germline/03_genotyping.sh \\
        --gvcf-args "${gvcf_args}" \\
        --outdir .
    """
}

// ── Process 4: Hard Filter + Annotate ────────────────────────────────────────
// Replaces VQSR (which requires large callsets) with hard filters on GATK INFO
// annotations (QD, MQ, FS, SOR, ReadPosRankSum, MappingQualityRankSum).
// Then annotates passing variants with the full annotation stack.
process HARD_FILTER_ANNOTATE {
    label 'gatk'

    publishDir "${params.outdir}/germline/filtered",  mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/germline/annotated", mode: 'copy', pattern: '*.tsv'

    input:
    path vcf
    path tbi

    output:
    path "hardfiltered.vcf.gz",     emit: vcf
    path "hardfiltered.vcf.gz.tbi", emit: tbi
    path "annotated.tsv",           emit: tsv

    script:
    """
    bash ${projectDir}/germline/04_hardfilter_annotate.sh \\
        --vcf ${vcf} \\
        --outdir .

    python ${projectDir}/annotation/annotate_variants.py \\
        --vcf hardfiltered.vcf.gz \\
        --type germline \\
        --output annotated.tsv
    """
}

// ── Germline workflow ─────────────────────────────────────────────────────────
workflow GERMLINE_WORKFLOW {

    samples = Channel.of(
        tuple('mother', file(cfg.paths.germline.bam_mother), file(cfg.paths.germline.bam_mother + ".bai")),
        tuple('father', file(cfg.paths.germline.bam_father), file(cfg.paths.germline.bam_father + ".bai")),
        tuple('son',    file(cfg.paths.germline.bam_son),    file(cfg.paths.germline.bam_son    + ".bai"))
    )

    MARK_DUPLICATES_GERMLINE(samples)

    HAPLOTYPECALLER_GVCF(MARK_DUPLICATES_GERMLINE.out.bam)

    gvcfs = HAPLOTYPECALLER_GVCF.out.map { sample, gvcf, tbi -> gvcf }.collect()
    tbis  = HAPLOTYPECALLER_GVCF.out.map { sample, gvcf, tbi -> tbi  }.collect()

    JOINT_GENOTYPING(gvcfs, tbis)

    HARD_FILTER_ANNOTATE(
        JOINT_GENOTYPING.out.vcf,
        JOINT_GENOTYPING.out.tbi
    )
}
