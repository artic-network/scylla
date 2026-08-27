#!/usr/bin/env python

import pyfastx
import argparse
from collections import defaultdict
import sys

from assignment import trim_read_id


def check_fastq(read_file, sample_size=1000):
    is_duplicates = True
    is_interleaved = False
    is_concat = False

    position = 0
    # Reference structures are only ever populated from the first
    # `sample_size` reads, so memory is bounded regardless of file size.
    # Reads beyond that window are only *probed* against this reference,
    # never added to it - any duplicate whose first occurrence falls
    # outside the window is missed, but for both the interleaved case
    # (pairs are adjacent) and the concatenated case (the whole first
    # half is duplicated verbatim in the second half) every duplicate's
    # first occurrence is guaranteed to fall within an early window.
    ref_names = set()
    ref_trimmed = set()
    ref_positions = {}
    ref_checks = {}

    differences = defaultdict(int)
    dup_positions = []
    early_exit = False

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record
        trimmed_name = trim_read_id(name)

        if position < sample_size:
            if name in ref_names or trimmed_name in ref_trimmed:
                differences[trimmed_name] = position - ref_positions[trimmed_name]
                dup_positions.append(position)
            if trimmed_name in ref_checks and ref_checks[trimmed_name] != seq:
                is_duplicates = False
            ref_positions[trimmed_name] = position
            ref_names.add(name)
            ref_trimmed.add(trimmed_name)
            if str(position + 1).startswith("1") or str(position + 1).startswith("5"):
                ref_checks[trimmed_name] = seq
        elif trimmed_name in ref_trimmed:
            differences[trimmed_name] = position - ref_positions[trimmed_name]
            dup_positions.append(position)
            if trimmed_name in ref_checks and ref_checks[trimmed_name] != seq:
                is_duplicates = False
            ref_positions[trimmed_name] = position

        position += 1

        # Once we're well past the reference window and every duplicate
        # found so far is an adjacent pair, this is confidently an
        # interleaved file - stop scanning early rather than reading the
        # rest of a potentially very large FASTQ just to confirm it.
        if (
            position >= 2 * sample_size
            and differences
            and set(differences.values()) == {1}
        ):
            early_exit = True
            break

    num_seqs = position

    # if no duplicated names or trimmed_names, then no need to do anything
    if len(differences) == 0:
        return 0

    difference_set = set([v for v in differences.values()])
    min_duplicate = min(dup_positions)
    if difference_set == {1}:
        # if all pairs are next to each other, have interleaved file
        is_interleaved = True
    else:
        # otherwise assume concatenated file
        is_concat = True

    sys.stderr.write(
        f"Found evidence of interleaving: {is_interleaved}, concatenation: {is_concat}, duplicates: {is_duplicates}.\nSplitting out reads\n"
    )
    out_prefix = read_file.split("/")[-1].split(".")[0]

    counts = defaultdict(int)
    if is_duplicates:
        sys.stderr.write(f"Position of first duplicate: {min_duplicate}\n")
        position = 0
        with open(f"{out_prefix}.fixed.fastq", "w") as r:
            for record in pyfastx.Fastq(read_file, build_index=False):
                name, seq, qual = record
                if position < min_duplicate:
                    r.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    counts["r"] += 1
                else:
                    trimmed_name = trim_read_id(name)
                    assert trimmed_name in differences
                position += 1
                if position >= min_duplicate:
                    break

    elif is_interleaved:
        with open(f"{out_prefix}.R1.fastq", "w") as r1, open(
            f"{out_prefix}.R2.fastq", "w"
        ) as r2:
            last = None
            for record in pyfastx.Fastq(read_file, build_index=False):
                name, seq, qual = record
                trimmed_name = trim_read_id(name)
                if last and trimmed_name == last:
                    r2.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    counts["r2"] += 1
                else:
                    r1.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    counts["r1"] += 1
                last = trimmed_name

    elif is_concat:
        position = 0
        with open(f"{out_prefix}.R1.fastq", "w") as r1, open(
            f"{out_prefix}.R2.fastq", "w"
        ) as r2:
            out_handle = r1
            key = "r1"
            for record in pyfastx.Fastq(read_file, build_index=False):
                name, seq, qual = record
                trimmed_name = trim_read_id(name)
                out_handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
                counts[key] += 1
                position += 1
                if position >= min_duplicate:
                    out_handle = r2
                    key = "r2"

    num_seqs_desc = f">= {num_seqs}" if early_exit else str(num_seqs)
    if is_duplicates:
        sys.stderr.write(
            f"Input {num_seqs_desc} sequences have resulted in out file with the following read counts: {out_prefix}.fixed.fastq : {counts['r']}\n"
        )
    else:
        sys.stderr.write(
            f"Input {num_seqs_desc} sequences have resulted in out files with the following read counts: {out_prefix}.R1.fastq : {counts['r1']}, {out_prefix}.R2.fastq : {counts['r2']}\n"
        )

    return 11


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Check a single FASTQ file to see if it contains paired reads, and split if it does."
        )
    )
    parser.add_argument("--fastq", help="Input FASTQ.")
    parser.add_argument(
        "--sample-size",
        dest="sample_size",
        type=int,
        default=1000,
        help=(
            "Number of leading reads to use to build the reference set for "
            "duplicate/interleave/concatenation detection (default: 1000)."
        ),
    )

    args = parser.parse_args()

    exit_code = check_fastq(args.fastq, sample_size=args.sample_size)
    sys.exit(exit_code)
