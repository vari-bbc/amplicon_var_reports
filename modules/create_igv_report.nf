process CREATE_IGV_REPORT {

    tag "${sample_id}:${amplicon}"

    label 'process_low'

    conda '/varidata/research/projects/bbc/research/KEMT_20260330_VBCS-1659/conda-1659'

    publishDir "${params.outdir}/${sample_id}/${amplicon}/igv_report", mode: "copy"

    input:
    tuple val(sample_id), val(amplicon), path(reference), path(bam), path(bai), path(vcf), path(tbi)

    output:
    tuple val(sample_id), val(amplicon), path("${sample_id}.${amplicon}.igv_report.html")

    script:
    def flanking = params.igv_report_flanking != null ? params.igv_report_flanking : 1000
    def standalone_arg = params.igv_report_standalone ? '--standalone' : ''
    """
    if [[ ! -f "${reference}.fai" ]]; then
      samtools faidx ${reference}
    fi

    create_report ${vcf} \\
      --fasta ${reference} \\
      --tracks ${vcf} ${bam} \\
      --flanking ${flanking} \\
      --title "${sample_id}:${amplicon} IGV Report" \\
      ${standalone_arg} \\
      --output ${sample_id}.${amplicon}.igv_report.html
    """
}
