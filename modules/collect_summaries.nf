process COLLECT_SUMMARIES {

    publishDir "${params.outdir}/summary", mode: "copy"

    input:
    path(classification_summaries)
    path(variant_summaries)

    output:
    path("all_amplicon_read_classification_summary.tsv")
    path("all_amplicon_variant_summary.tsv")

    script:
    """
    head -n 1 ${classification_summaries[0]} > all_amplicon_read_classification_summary.tsv
    for f in ${classification_summaries}; do
      tail -n +2 \$f >> all_amplicon_read_classification_summary.tsv
    done

    head -n 1 ${variant_summaries[0]} > all_amplicon_variant_summary.tsv
    for f in ${variant_summaries}; do
      tail -n +2 \$f >> all_amplicon_variant_summary.tsv
    done
    """
}