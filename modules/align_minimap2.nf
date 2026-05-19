process ALIGN_MINIMAP2 {

    tag "${sample_id}:${amplicon}"

    module 'bbc2/minimap2/minimap2-2.30'
    module 'bbc2/samtools/samtools-1.21'

    publishDir "${params.outdir}/${sample_id}/${amplicon}/bam", mode: "copy"

    input:
    tuple val(sample_id), val(amplicon), path(fastq), path(reference)

    output:
    tuple val(sample_id), val(amplicon), path(reference), path("${sample_id}.${amplicon}.bam"), path("${sample_id}.${amplicon}.bam.bai")

    script:
    """
    minimap2 -t ${params.threads} -ax map-ont ${reference} ${fastq} | \\
      samtools sort -@ ${params.threads} -o ${sample_id}.${amplicon}.bam

    samtools index ${sample_id}.${amplicon}.bam
    """
}
