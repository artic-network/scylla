#!/usr/bin/env python

import sys
import argparse
import pyfastx
import json


def main():
    # Parse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s",
        "-s1",
        "-1",
        dest="reads1",
        required=True,
        help="FASTA/FASTQ File containing the raw reads.",
    )
    parser.add_argument(
        "-s2",
        "-2",
        dest="reads2",
        default="",
        help="2nd FASTA/FASTQ File containing the raw reads (paired).",
    )

    parser.set_defaults(append=False)
    args = parser.parse_args()

    fq = pyfastx.Fastq(args.reads1)
    total_length = fq.size

    if (args.reads2):
        fq = pyfastx.Fastq(args.reads2)
        total_length += fq.size

    with open(f"total_length.json", "w") as f:
            json.dump({"total_len": total_length}, f)

    sys.exit(0)


if __name__ == "__main__":
    main()
