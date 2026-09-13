#!/usr/bin/env python
"""
No dependency Python script for joining two paired end FASTQ files.
Supports concatenating reads with a separator ("NNNNN") or interleaving reads via the
--interleave option. Auto-detects gzip'd files, offers header checking via a --strict flag,
and supports output to STDOUT a gzip'd FASTQ or an uncompressed FASTQ (--uncompressed flag).
"""
import argparse
import sys
import gzip


def enforce_headers(f1_header, f2_header):
    if f1_header[0] != "@" or f2_header[0] != "@":
        sys.exit(5)

    if f1_header.split()[0].endswith("/1") and f2_header.split()[0].endswith("/2"):
        if f1_header.split()[0][:-2] != f2_header.split()[0][:-2]:
            sys.exit(8)

    elif f1_header.strip().split(" ")[0] != f2_header.strip().split(" ")[0]:
        sys.exit(8)


FLUSH_BYTES = 1 << 20


def _close_out(f_out):
    if f_out is not sys.stdout:
        f_out.close()


def join_w_separator(f1, f2, f_out, strict=False, sep="|"):
    ix = 0

    # Use 'G', which is in most Phred scores.
    # http://en.wikipedia.org/wiki/FASTQ_format
    phred_sep = "|" * len(sep)

    buf = []
    buf_len = 0

    while True:
        f1_line = f1.readline()
        f2_line = f2.readline()

        if f1_line == "":
            if buf:
                f_out.write("".join(buf))
            f1.close()
            f2.close()
            _close_out(f_out)
            if f2_line != "":
                raise Exception("Input FASTQ files do not match in length.")
            return

        if ix % 4 == 0:  # Header
            if strict:  # Fail if they don't match up to the first whitespace
                enforce_headers(f1_line, f2_line)

            # Write the header out
            line = f1_line

        elif (ix - 1) % 4 == 0:  # Sequence
            line = f1_line.strip() + sep + f2_line
        elif (ix - 2) % 4 == 0:  # Separator
            line = "+\n"
        else:
            line = f1_line.strip() + phred_sep + f2_line

        buf.append(line)
        buf_len += len(line)
        if buf_len >= FLUSH_BYTES:
            f_out.write("".join(buf))
            buf = []
            buf_len = 0

        ix += 1


def join_interleaved(f1, f2, f_out, strict=False):
    ix = 0
    f1_lines = []
    f2_lines = []
    buf = []
    buf_len = 0

    def flush_record(f1_lines, f2_lines):
        nonlocal buf, buf_len
        if f1_lines and f2_lines:
            assert len(f1_lines) == 4
            assert len(f2_lines) == 4
            for line in f1_lines:
                buf.append(line)
                buf_len += len(line)
            for line in f2_lines:
                buf.append(line)
                buf_len += len(line)

            f1_lines = []
            f2_lines = []

        return f1_lines, f2_lines

    while True:
        f1_line = f1.readline()
        f2_line = f2.readline()

        if f1_line == "":
            if f2_line != "":
                raise Exception("Input FASTQ files do not match in length.")
            break

        if ix % 4 == 0:  # Header #1
            if strict:
                enforce_headers(f1_line, f2_line)

            f1_lines, f2_lines = flush_record(f1_lines, f2_lines)
            if buf_len >= FLUSH_BYTES:
                f_out.write("".join(buf))
                buf = []
                buf_len = 0

        # Fill buffer up to 4 lines
        f1_lines.append(f1_line)
        f2_lines.append(f2_line)
        ix += 1

    flush_record(f1_lines, f2_lines)
    if buf:
        f_out.write("".join(buf))
    f1.close()
    f2.close()
    _close_out(f_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Join two FASTQ paired end read files. Defaults to interleaving reads. "
            "Note: Requires input files to be sorted."
        )
    )
    parser.add_argument("fastq1", help="First input FASTQ.")
    parser.add_argument("fastq2", help="Second input FASTQ.")
    parser.add_argument(
        "output_fastq",
        nargs="?",
        help=("Output FASTQ file name (optional, " "streams to STDOUT if missing."),
    )
    parser.add_argument(
        "--sep",
        default="|",
        help=("Optional separator to override default (|)."),
    )
    parser.add_argument(
        "--no-interleave",
        action="store_true",
        help="Concatenate the two reads using a separate (--sep) rather than "
        "interleaving them (default).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Enforce that the headers of the input "
            "FASTQ files match until the first whitespace "
            "delimiter (e.g., a space)."
        ),
    )
    parser.add_argument("--gzip", action="store_true", help=("Gzip output_fastq."))
    args = parser.parse_args()

    if ".gz" in args.fastq1 or ".gzip" in args.fastq1:
        f1 = gzip.open(args.fastq1, mode="rt")
    else:
        f1 = open(args.fastq1, mode="rt", buffering=FLUSH_BYTES)
    if ".gz" in args.fastq2 or ".gzip" in args.fastq2:
        f2 = gzip.open(args.fastq2, mode="rt")
    else:
        f2 = open(args.fastq2, mode="rt", buffering=FLUSH_BYTES)

    if args.output_fastq is None:
        f_out = sys.stdout
    elif args.gzip:
        f_out = gzip.open(args.output_fastq, mode="wt")
    else:
        f_out = open(args.output_fastq, mode="w", buffering=FLUSH_BYTES)

    if not args.no_interleave:
        join_interleaved(f1, f2, f_out, strict=args.strict)
    else:
        join_w_separator(f1, f2, f_out, strict=args.strict, sep=args.sep)
