process CALL_VARIANTS {

    tag "${sample_id}:${amplicon}"

    module 'bbc2/samtools/samtools-1.21'
    module 'bbc2/bcftools/bcftools-1.23'

    publishDir "${params.outdir}/${sample_id}/${amplicon}/vcf", mode: "copy"

    input:
    tuple val(sample_id), val(amplicon), path(reference), path(bam), path(bai)

    output:
    tuple val(sample_id), val(amplicon), path(reference), path("${sample_id}.${amplicon}.vcf.gz"), path("${sample_id}.${amplicon}.vcf.gz.tbi")

    script:
    """
    bcftools mpileup \\
      -Ou \\
      -f ${reference} \\
      ${bam} | \\
    bcftools call \\
      -mv \\
      -Oz \\
      -o ${sample_id}.${amplicon}.vcf.gz

    bcftools index -t ${sample_id}.${amplicon}.vcf.gz
    """
}