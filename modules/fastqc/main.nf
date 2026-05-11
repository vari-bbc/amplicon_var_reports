process FASTQC {
    tag "$fq_pref"
    label 'process_medium'
    module 'bbc2/fastqc/fastqc-0.12.1'
  input:
    tuple val(fq_pref), path("${fq_pref}.fastq.gz")

  output:
    path "${fq_pref}_fastqc.html", emit: html
    path "${fq_pref}_fastqc.zip", emit: zip

  script:
    """
    fastqc "${fq_pref}.fastq.gz"
    """
}
