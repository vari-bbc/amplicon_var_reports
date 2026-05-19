process CLASSIFY_READS {

    tag "${sample_id}:${amplicon}"

    conda '/varidata/research/projects/bbc/research/KEMT_20260330_VBCS-1659/conda-1659'

    publishDir "${params.outdir}/${sample_id}/${amplicon}/classification", mode: "copy"

    input:
    tuple val(sample_id), val(amplicon), path(reference), path(bam), path(bai)

    output:
    path("${sample_id}.${amplicon}.read_classification.tsv")
    path("${sample_id}.${amplicon}.read_classification.summary.tsv")

    script:
    """
    python ${workflow.projectDir}/scripts/classify_reads.py \\
      --sample-id ${sample_id} \\
      --amplicon ${amplicon} \\
      --reference ${reference} \\
      --bam ${bam} \\
      --min-mapq ${params.min_mapq} \\
      --min-ref-cov-frac ${params.min_ref_cov_frac} \\
      --min-identity ${params.min_identity} \\
      --max-indel-bp ${params.max_indel_bp} \\
      --max-softclip-frac ${params.max_softclip_frac} \\
      --per-read-out ${sample_id}.${amplicon}.read_classification.tsv \\
      --summary-out ${sample_id}.${amplicon}.read_classification.summary.tsv
    """
}
