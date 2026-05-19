#!/usr/bin/env python3

import argparse
import csv
from collections import Counter

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify Nanopore amplicon reads as expected/correct or mutated based on alignment to expected reference."
    )

    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--amplicon", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--bam", required=True)

    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-ref-cov-frac", type=float, default=0.90)
    parser.add_argument("--min-identity", type=float, default=0.98)
    parser.add_argument("--max-indel-bp", type=int, default=10)
    parser.add_argument("--max-softclip-frac", type=float, default=0.20)

    parser.add_argument("--per-read-out", required=True)
    parser.add_argument("--summary-out", required=True)

    return parser.parse_args()


def cigar_counts(read):
    matches = 0
    insertions = 0
    deletions = 0
    softclips = 0
    hardclips = 0
    skipped = 0

    if read.cigartuples is None:
        return {
            "matches": 0,
            "insertions": 0,
            "deletions": 0,
            "softclips": 0,
            "hardclips": 0,
            "skipped": 0,
        }

    for op, length in read.cigartuples:
        if op in {0, 7, 8}:      # M, =, X
            matches += length
        elif op == 1:            # I
            insertions += length
        elif op == 2:            # D
            deletions += length
        elif op == 4:            # S
            softclips += length
        elif op == 5:            # H
            hardclips += length
        elif op == 3:            # N
            skipped += length

    return {
        "matches": matches,
        "insertions": insertions,
        "deletions": deletions,
        "softclips": softclips,
        "hardclips": hardclips,
        "skipped": skipped,
    }


def count_mismatches_against_reference(read, ref_seqs):
    """
    Count mismatches by comparing aligned query bases to reference bases.

    This avoids relying on NM/MD tags.
    """
    if read.is_unmapped:
        return 0, 0

    ref_name = read.reference_name
    mismatches = 0
    aligned_bases = 0

    pairs = read.get_aligned_pairs(matches_only=False)

    for query_pos, ref_pos in pairs:
        if query_pos is None or ref_pos is None:
            continue

        query_base = read.query_sequence[query_pos].upper()
        ref_base = ref_seqs[ref_name][ref_pos].upper()

        if query_base not in {"A", "C", "G", "T"}:
            continue
        if ref_base not in {"A", "C", "G", "T"}:
            continue

        aligned_bases += 1

        if query_base != ref_base:
            mismatches += 1

    return mismatches, aligned_bases


def classify_read(
    read,
    ref_seqs,
    ref_lengths,
    min_mapq,
    min_ref_cov_frac,
    min_identity,
    max_indel_bp,
    max_softclip_frac,
):
    if read.is_unmapped:
        return "unmapped", {}

    if read.is_secondary or read.is_supplementary:
        return "non_primary_alignment", {}

    if read.mapping_quality < min_mapq:
        return "low_mapq", {}

    ref_name = read.reference_name
    ref_len = ref_lengths[ref_name]

    ref_covered = read.reference_end - read.reference_start
    ref_cov_frac = ref_covered / ref_len if ref_len > 0 else 0

    counts = cigar_counts(read)
    mismatches, aligned_bases = count_mismatches_against_reference(read, ref_seqs)

    insertion_bp = counts["insertions"]
    deletion_bp = counts["deletions"]
    indel_bp = insertion_bp + deletion_bp
    softclip_bp = counts["softclips"]

    query_len = read.query_length or 0
    softclip_frac = softclip_bp / query_len if query_len > 0 else 1

    total_error_bp = mismatches + insertion_bp + deletion_bp
    identity = 1 - (total_error_bp / aligned_bases) if aligned_bases > 0 else 0

    metrics = {
        "ref_name": ref_name,
        "ref_start_0based": read.reference_start,
        "ref_end_0based": read.reference_end,
        "ref_covered_bp": ref_covered,
        "ref_length_bp": ref_len,
        "ref_cov_frac": ref_cov_frac,
        "mapq": read.mapping_quality,
        "query_length": query_len,
        "aligned_bases": aligned_bases,
        "mismatches": mismatches,
        "insertion_bp": insertion_bp,
        "deletion_bp": deletion_bp,
        "indel_bp": indel_bp,
        "softclip_bp": softclip_bp,
        "softclip_frac": softclip_frac,
        "identity": identity,
    }

    if ref_cov_frac < min_ref_cov_frac:
        return "partial_alignment", metrics

    if softclip_frac > max_softclip_frac:
        return "large_softclip_or_extra_sequence", metrics

    if indel_bp > max_indel_bp:
        return "indel_mutated", metrics

    if identity < min_identity:
        return "substitution_or_error_rich", metrics

    return "correct_expected", metrics


def main():
    args = parse_args()

    fasta = pysam.FastaFile(args.reference)
    bam = pysam.AlignmentFile(args.bam, "rb")

    ref_lengths = {
        ref: fasta.get_reference_length(ref)
        for ref in fasta.references
    }

    ref_seqs = {
        ref: fasta.fetch(ref).upper()
        for ref in fasta.references
    }

    fieldnames = [
        "sample_id",
        "amplicon",
        "read_id",
        "classification",
        "ref_name",
        "ref_start_0based",
        "ref_end_0based",
        "ref_covered_bp",
        "ref_length_bp",
        "ref_cov_frac",
        "mapq",
        "query_length",
        "aligned_bases",
        "mismatches",
        "insertion_bp",
        "deletion_bp",
        "indel_bp",
        "softclip_bp",
        "softclip_frac",
        "identity",
    ]

    counts = Counter()
    total_primary_reads_seen = 0

    with open(args.per_read_out, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue

            total_primary_reads_seen += 1

            classification, metrics = classify_read(
                read=read,
                ref_seqs=ref_seqs,
                ref_lengths=ref_lengths,
                min_mapq=args.min_mapq,
                min_ref_cov_frac=args.min_ref_cov_frac,
                min_identity=args.min_identity,
                max_indel_bp=args.max_indel_bp,
                max_softclip_frac=args.max_softclip_frac,
            )

            counts[classification] += 1

            row = {
                "sample_id": args.sample_id,
                "amplicon": args.amplicon,
                "read_id": read.query_name,
                "classification": classification,
            }

            for key in fieldnames:
                if key not in row:
                    row[key] = metrics.get(key, "")

            writer.writerow(row)

    with open(args.summary_out, "w", newline="") as summary_handle:
        fieldnames_summary = [
            "sample_id",
            "amplicon",
            "classification",
            "read_count",
            "total_reads",
            "percent_reads",
        ]

        writer = csv.DictWriter(summary_handle, fieldnames=fieldnames_summary, delimiter="\t")
        writer.writeheader()

        for classification, count in sorted(counts.items()):
            pct = 100 * count / total_primary_reads_seen if total_primary_reads_seen else 0

            writer.writerow({
                "sample_id": args.sample_id,
                "amplicon": args.amplicon,
                "classification": classification,
                "read_count": count,
                "total_reads": total_primary_reads_seen,
                "percent_reads": f"{pct:.3f}",
            })

    bam.close()
    fasta.close()


if __name__ == "__main__":
    main()