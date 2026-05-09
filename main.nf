#!/usr/bin/env nextflow
// =============================================================================
// main.nf  —  Top-level Nextflow DSL2 workflow
//
// Orchestrates the somatic and/or germline variant calling modules.
// Each module is a self-contained Nextflow module file in modules/.
//
// Usage:
//   nextflow run main.nf -profile local --module somatic
//   nextflow run main.nf -profile local --module germline
//   nextflow run main.nf -profile local --module all
//
// What Nextflow adds over plain bash scripts:
//   - Automatic dependency tracking: only re-runs steps whose inputs changed
//   - Parallel execution of independent steps (e.g. 3 HaplotypeCaller jobs)
//   - Structured output routing via publishDir
//   - HTML execution report + timeline for every run
//   - Retry logic on transient failures
// =============================================================================

nextflow.enable.dsl = 2

include { SOMATIC_WORKFLOW } from './modules/somatic'
include { GERMLINE_WORKFLOW } from './modules/germline'

// ── Parameter validation ──────────────────────────────────────────────────────
def valid_modules = ['somatic', 'germline', 'all']
if (!valid_modules.contains(params.module)) {
    error "Unknown --module '${params.module}'. Choose from: ${valid_modules.join(', ')}"
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    log.info """
    ╔═══════════════════════════════════════════════════════╗
    ║  Variant Calling Pipeline                             ║
    ║  Module  : ${params.module.padRight(43)}║
    ║  Config  : ${params.config.padRight(43)}║
    ║  Output  : ${params.outdir.padRight(43)}║
    ╚═══════════════════════════════════════════════════════╝
    """.stripIndent()

    if (params.module == 'somatic' || params.module == 'all') {
        SOMATIC_WORKFLOW()
    }

    if (params.module == 'germline' || params.module == 'all') {
        GERMLINE_WORKFLOW()
    }
}

// ── Completion handler ────────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    log.info """
    Pipeline ${status}
    Duration  : ${workflow.duration}
    Results   : ${params.outdir}/
    Report    : ${params.tracedir}/report.html
    """.stripIndent()
}
