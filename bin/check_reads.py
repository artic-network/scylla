#!/usr/bin/env python

import pyfastx
import argparse
from collections import defaultdict
import sys

from assignment import trim_read_id


def check_fastq(
    read_file, sample_size=1000, max_illumina_read_length=1000, full_scan=True
):
    """
    Detect whether a single FASTQ actually holds paired reads - interleaved,
    concatenated, or straightforwardly duplicated - and split it if so.

    Args:
        read_file (str): Path to the FASTQ to check.
        sample_size (int): Number of leading reads used to build the reference
            set of read names. Ignored for name tracking when full_scan is set,
            where every read joins the reference; still bounds the much more
            expensive sequence snapshots either way.
        max_illumina_read_length (int): Read length above which the file cannot
            be Illumina paired-end data. Only consulted when full_scan is False.
        full_scan (bool): Examine the whole file - the default. Every read is
            read, and every read joins the reference set, so a duplicate is
            found wherever its first occurrence falls. Clearing this enables
            three shortcuts, each of which trades detection for throughput:

              * a read longer than max_illumina_read_length proves the file is
                not Illumina paired-end and abandons the check immediately. That
                misses a long-read file concatenated onto itself, because the
                duplicate evidence only appears once the scan reaches the
                original file's length - by which point a long read has almost
                always already been seen;
              * once past the reference window with every duplicate so far an
                adjacent pair, the file is confidently interleaved and scanning
                stops. That misses any non-adjacent duplicate later in the file,
                and makes the reported sequence count a lower bound;
              * only the leading `sample_size` reads join the reference, so any
                duplicate whose first occurrence falls beyond the window is
                invisible. This is what keeps memory flat regardless of input
                size, and it costs nothing for the two cases this check targets,
                since interleaved pairs are adjacent and a concatenated file
                repeats its whole first half.

    The cost of the default is that the reference grows with the file: one dict
    entry per distinct trimmed read name. Sequence snapshots stay bounded to the
    leading window in both modes, since holding a sequence per read would
    dominate that on long-read input.

    Returns:
        int: 0 if nothing needed fixing, 11 if the file was split.
    """
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
    # Trimmed read name -> the position it was last seen at. This is the whole
    # reference: a previously kept set of full names and a separate set of
    # trimmed names were both redundant with it. trim_read_id is deterministic,
    # so seeing a full name again implies seeing its trimmed form again, and the
    # set of trimmed names is exactly this dict's key set. Collapsing the three
    # matters because under a full scan this grows with the file.
    ref_positions = {}
    # Trimmed read name -> sequence, for a sparse sample of early reads. Used
    # only to tell real mates (different sequences) from a verbatim duplicate.
    ref_checks = {}

    differences = defaultdict(int)
    dup_positions = []
    early_exit = False

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record

        # Illumina chemistry has a hard per-read length ceiling well under
        # 1kb (even long-read kits top out around 600bp), so a single read
        # longer than this proves the file cannot be interleaved or
        # concatenated Illumina paired-end data - the entire pairing check
        # is moot and we can stop scanning immediately, which matters a lot
        # for very large long-read (e.g. ONT) FASTQs. Checked on every read
        # for the whole scan (not just the reference window) since it's a
        # free len() on data already parsed, and a degraded/fragmented run
        # could have a run of short reads before the long ones appear.
        #
        # Gated on `not differences` so it never fires mid-detection once
        # duplicate evidence has actually been seen. In practice this rarely
        # protects the "whole ONT file accidentally duplicated onto itself"
        # case though: that duplicate evidence only appears once the scan
        # reaches the original file's length, and a long-read file almost
        # always has a read over the threshold well before that - so this
        # optimisation trades away that (rarer, non-Illumina-specific)
        # detection in favour of the much more common speed win. That trade is
        # exactly what full_scan declines to make.
        if not full_scan and not differences and len(seq) > max_illumina_read_length:
            sys.stderr.write(
                f"Read {name} is {len(seq)}bp (> {max_illumina_read_length}bp): "
                "file cannot be Illumina paired-end data, skipping pairing check\n"
            )
            return 0

        trimmed_name = trim_read_id(name)

        # Under a full scan every read joins the reference, so a duplicate is
        # caught wherever in the file its first occurrence falls. Otherwise only
        # the leading `sample_size` reads do: later reads are still probed
        # against that reference but never added to it, which bounds memory at
        # the cost of missing any duplicate first seen past the window.
        if full_scan or position < sample_size:
            if trimmed_name in ref_positions:
                differences[trimmed_name] = position - ref_positions[trimmed_name]
                dup_positions.append(position)
            if trimmed_name in ref_checks and ref_checks[trimmed_name] != seq:
                is_duplicates = False
            ref_positions[trimmed_name] = position
            # Sequence snapshots stay bounded to the leading window even under a
            # full scan. They only need to catch mates whose sequences differ,
            # and keeping a sequence per read would dominate memory on long-read
            # input far more than the names do.
            if position < sample_size and (
                str(position + 1).startswith("1") or str(position + 1).startswith("5")
            ):
                ref_checks[trimmed_name] = seq
        elif trimmed_name in ref_positions:
            differences[trimmed_name] = position - ref_positions[trimmed_name]
            dup_positions.append(position)
            if trimmed_name in ref_checks and ref_checks[trimmed_name] != seq:
                is_duplicates = False
            ref_positions[trimmed_name] = position

        position += 1

        # Once we're well past the reference window and every duplicate
        # found so far is an adjacent pair, this is confidently an
        # interleaved file - stop scanning early rather than reading the
        # rest of a potentially very large FASTQ just to confirm it. Under
        # full_scan we keep reading, so a non-adjacent duplicate further in
        # still reclassifies the file as concatenated rather than interleaved.
        if (
            not full_scan
            and position >= 2 * sample_size
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
    parser.add_argument(
        "--max-illumina-read-length",
        dest="max_illumina_read_length",
        type=int,
        default=1000,
        help=(
            "A single read longer than this (bp) proves the file cannot be "
            "Illumina paired-end data, so the pairing check is skipped "
            "immediately (default: 1000)."
        ),
    )

    parser.add_argument(
        "--no-full-scan",
        dest="no_full_scan",
        action="store_true",
        default=False,
        help=(
            "Stop reading as soon as the answer looks settled, instead of "
            "scanning every record. Faster on large inputs, at the cost of "
            "missing a long-read file concatenated onto itself and any "
            "non-adjacent duplicate that appears late in an otherwise "
            "interleaved file."
        ),
    )

    args = parser.parse_args()

    exit_code = check_fastq(
        args.fastq,
        sample_size=args.sample_size,
        max_illumina_read_length=args.max_illumina_read_length,
        full_scan=not args.no_full_scan,
    )
    sys.exit(exit_code)
