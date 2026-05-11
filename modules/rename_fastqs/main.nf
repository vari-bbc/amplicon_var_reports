process RENAME_FASTQS {
    tag "$sample"
    label 'process_medium'
    input:
    tuple val(sample), path(fastq), val(suffix)

    output:
    path "${sample}_${suffix}.fastq.gz"

    script:
    if (fastq.size() == 1) {
        """
        ln -s ${fastq[0]} "${sample}_${suffix}.fastq.gz"

        """
    } else {
        """
        cat ${fastq.join(' ')} > "${sample}_${suffix}.fastq.gz"

        """
    }
}
