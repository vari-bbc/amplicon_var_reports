# VBCS-1659

The pipeline writes per-sample, per-amplicon IGV reports to:

`results/<sample_id>/<amplicon>/igv_report/<sample_id>.<amplicon>.igv_report.html`

The IGV report module uses the called VCF as the site table and includes both the VCF and BAM as IGV tracks. Report padding can be changed with `--igv_report_flanking`. Set `--igv_report_standalone true` to ask `igv-reports` to embed JavaScript in the HTML report.

## Example

```shell
nextflow run main.nf -profile slurm --samplesheet samplesheet.tsv --threads 8 
```