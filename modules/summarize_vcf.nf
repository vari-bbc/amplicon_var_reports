process SUMMARIZE_VCF {

    tag "${sample_id}:${amplicon}"

    conda '/varidata/research/projects/bbc/research/KEMT_20260330_VBCS-1659/conda-1659'

    publishDir "${params.outdir}/${sample_id}/${amplicon}/vcf_summary", mode: "copy"

    input:
    tuple val(sample_id), val(amplicon), path(reference), path(vcf), path(tbi)

    output:
    path("${sample_id}.${amplicon}.variant_summary.tsv")

    script:
    """
    python ${workflow.projectDir}/scripts/summarize_vcf.py \\
      --sample-id ${sample_id} \\
      --amplicon ${amplicon} \\
      --vcf ${vcf} \\
      --out ${sample_id}.${amplicon}.variant_summary.tsv
    """
}