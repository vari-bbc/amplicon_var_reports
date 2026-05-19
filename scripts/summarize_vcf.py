#!/usr/bin/env python3

import argparse
import csv
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize variants from a VCF."
    )

    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--amplicon", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)

    return parser.parse_args()


def infer_variant_type(ref, alts):
    if not alts:
        return "no_alt"

    alt = alts[0]

    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if len(ref) < len(alt):
        return "INS"

    if len(ref) > len(alt):
        return "DEL"

    return "MNV_or_complex"


def get_info_value(record, key):
    if key not in record.info:
        return ""

    value = record.info[key]
    if value is None:
        return ""

    if isinstance(value, tuple):
        return ",".join(str(x) for x in value)

    return str(value)


def main():
    args = parse_args()

    vcf = pysam.VariantFile(args.vcf)

    fieldnames = [
        "sample_id",
        "amplicon",
        "chrom",
        "pos_1based",
        "ref",
        "alt",
        "variant_type",
        "qual",
        "filter",
        "depth_info_DP",
        "allele_depth_info_AD",
        "allele_freq_info_AF",
    ]

    with open(args.out, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for record in vcf.fetch():
            alts = record.alts or []

            writer.writerow({
                "sample_id": args.sample_id,
                "amplicon": args.amplicon,
                "chrom": record.chrom,
                "pos_1based": record.pos,
                "ref": record.ref,
                "alt": ",".join(alts),
                "variant_type": infer_variant_type(record.ref, alts),
                "qual": record.qual if record.qual is not None else "",
                "filter": ";".join(record.filter.keys()) if record.filter.keys() else "PASS",
                "depth_info_DP": get_info_value(record, "DP"),
                "allele_depth_info_AD": get_info_value(record, "AD"),
                "allele_freq_info_AF": get_info_value(record, "AF"),
            })

    vcf.close()


if __name__ == "__main__":
    main()