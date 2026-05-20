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

    if [[ "${vcf}" == *.gz ]]; then
      has_variants=\$(zgrep -v '^#' "${vcf}" | head -n 1 || true)
    else
      has_variants=\$(grep -v '^#' "${vcf}" | head -n 1 || true)
    fi

    if [[ -z "\${has_variants}" ]]; then
      echo "No variant records found in ${vcf}; using dummy site table to preserve coverage in IGV report"
      ref_name=\$(cut -f1 "${reference}.fai" | head -n 1)
      ref_len=\$(cut -f2 "${reference}.fai" | head -n 1)
      dummy_sites="${sample_id}.${amplicon}.dummy_sites.bed"
      printf "%s\t0\t%s\n" "\${ref_name}" "\${ref_len}" > "\${dummy_sites}"
      create_report "\${dummy_sites}" \
        --fasta ${reference} \
        --tracks ${bam} \
        --flanking ${flanking} \
        --title "${sample_id}:${amplicon} IGV Report" \
        ${standalone_arg} \
        --output ${sample_id}.${amplicon}.igv_report.html
    else
      create_report ${vcf} \
        --fasta ${reference} \
        --tracks ${vcf} ${bam} \
        --flanking ${flanking} \
        --title "${sample_id}:${amplicon} IGV Report" \
        ${standalone_arg} \
        --output ${sample_id}.${amplicon}.igv_report.html
    fi
    """
}
