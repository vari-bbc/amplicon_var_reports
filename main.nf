#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

params.samplesheet = params.samplesheet ?: "samplesheet.tsv"
params.outdir      = params.outdir ?: "results"


include { ALIGN_MINIMAP2 } from './modules/align_minimap2'
include { CALL_VARIANTS } from './modules/call_variants'
include { CLASSIFY_READS} from './modules/classify_reads'
include { SUMMARIZE_VCF} from './modules/summarize_vcf'
include { COLLECT_SUMMARIES} from './modules/collect_summaries'

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    sample_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id,
                row.amplicon,
                file(row.fastq),
                file(row.ref_fasta)
            )
        }
    mapped_ch = ALIGN_MINIMAP2(sample_ch)

    vcf_ch = CALL_VARIANTS(mapped_ch)

    (classification_ch, classification_summary_ch) = CLASSIFY_READS(mapped_ch)

    vcf_summary_ch = SUMMARIZE_VCF(vcf_ch)

    COLLECT_SUMMARIES(classification_ch, vcf_summary_ch)

    classification_summary_ch

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
