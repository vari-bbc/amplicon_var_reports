#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from './modules/fastqc'
include { RENAME_FASTQS } from './modules/rename_fastqs'
include { MULTIQC } from './modules/multiqc'


nextflow.preview.output = true

// Count letters for each word
process count_letters {
  input:
    val word

  output:
    path "${word}.txt"

  script:
    """
    printf "${word}" | wc -c > ${word}.txt
    """
}

// Sum all counts into one file
process sum_letter_counts {
  input:
    path(count_files)

  output:
    path "summed.txt"

  script:
    """
    cat ${count_files} | awk '{sum += \$1} END {print "Total:", sum}' > summed.txt
    """
}

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:

    Channel.fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> 
        tuple(row.sample, file("raw_data/" + row.fq1), file("raw_data/" + row.fq2))
    }
    .groupTuple()
        .multiMap { row -> 
        R1: [sample:row[0], fq:row[1], suff:"R1"]
        R2: [sample:row[0], fq:row[2], suff:"R2"]
    }
    .set { grouped_samples }

    renamed_fastqs = RENAME_FASTQS(grouped_samples.R1.concat(grouped_samples.R2))

    fastqc = FASTQC(RENAME_FASTQS
    .out
    .map { elem -> 
        def pref = elem.getBaseName().replace(".fastq", "")
        [pref:pref, fq:elem]
    }
    )

    multiqc = MULTIQC(
        FASTQC.out.html.
        map { elem -> elem.getParent() }
        .collect()
    )

    words_ch = Channel
              .fromPath('words.csv')
              .splitCsv(header: true)
              .map { item -> item['word'] }

    count_ch = count_letters(words_ch)
    sum_ch = sum_letter_counts(count_letters.out.collect())

    workflow.onComplete = {
        log.info(
            workflow.success
                ? "\nDone!\n"
                : "Oops .. something went wrong"
        )
    }

    publish:
    count_out = count_ch
    sum_out = sum_ch
    renamed_fastqs_out = renamed_fastqs
    fastqc_out = fastqc.html // "Cannot assign a multi-channel to a workflow output: fastqc_out"
    multiqc_out = multiqc

}

output {
    count_out {
        path 'step1'
    }
    sum_out {
        path 'step2'
    }
    renamed_fastqs_out {
        path 'renamed_fastqs'
    }
    fastqc_out {
        path 'fastqc'
    }
    multiqc_out {
        path 'multiqc'
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
