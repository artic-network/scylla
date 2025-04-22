#!/usr/bin/env python

import pyfastx
import argparse
from collections import defaultdict
import sys

from assignment import trim_read_id


def check_fastq(read_file):
    is_duplicates = True
    is_interleaved = False
    is_concat = False

    position = 0
    seen_names = set()
    seen_trimmed = set()
    positions = defaultdict(int)
    differences = defaultdict(int)
    checks = {}

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record
        trimmed_name = trim_read_id(name)

        if name in seen_names:
            assert trimmed_name not in differences
            differences[trimmed_name] = position - positions[trimmed_name]
        elif trimmed_name in seen_trimmed:
            differences[trimmed_name] = position - positions[trimmed_name]
        if trimmed_name in checks and checks[trimmed_name] != seq:
            is_duplicates = False
        positions[trimmed_name] = position
        position += 1
        seen_names.add(name)
        seen_trimmed.add(trimmed_name)
        if str(position).startswith("1") or str(position).startswith("5"):
            checks[trimmed_name] = seq

    num_seqs = position

    # if no duplicated names or trimmed_names, then no need to do anything
    if len(differences) == 0:
        return 0

    difference_set = set([v for v in differences.values()])
    position_set = set([positions[k] for k in differences])
    min_duplicate = min(position_set)
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

    if is_duplicates:
        sys.stderr.write(
            f"Input {num_seqs} sequences have resulted in out file with the following read counts: {out_prefix}.fixed.fastq : {counts['r']}\n"
        )
    else:
        sys.stderr.write(
            f"Input {num_seqs} sequences have resulted in out files with the following read counts: {out_prefix}.R1.fastq : {counts['r1']}, {out_prefix}.R2.fastq : {counts['r2']}\n"
        )

    return 11


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Check a single FASTQ file to see if it contains paired reads, and split if it does."
        )
    )
    parser.add_argument("--fastq", help="Input FASTQ.")

    args = parser.parse_args()

    exit_code = check_fastq(args.fastq)
    sys.exit(exit_code)
