#!/usr/bin/env python

#################
#
#   NOTE that the read counts/percentages are taken from the bracken reestimated file, but the reads themselves
#    are extracted based on kraken classifications because bracken does not provide classifications at the read level
#
###############

import sys
import os
import gzip
from Bio import SeqIO
import argparse
import json
from datetime import datetime
from collections import defaultdict

def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i != "":
            break
        depth += 1
    depth = int(depth / 2)
    return depth

def save_file(name, lines):
    if len(lines) == 0:
        return
    outfile = ".".join([name, "kreport_split.txt"])
    with open(outfile, "w") as f:
        for line in lines:
            f.write(line)

def parse_report_file(report_file, split_strings, split_rank):
    depth_dict = {}
    lines = {"remainder":[]}
    key = "remainder"
    hierarchy = []
    max_depth = -1

    # parses a kraken or bracken file
    with open(report_file, "r") as f:
        for line in f:
            if line.startswith("% of Seqs"):
                continue
            try:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")
            except:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    ignore1,
                    ignore2,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")
            depth = parse_depth(name)
            hierarchy = hierarchy[: depth]
            name = name.strip()
            add_hierarchy = False

            while depth <= max_depth:
                add_hierarchy = True
                del depth_dict[max_depth]
                if len(depth_dict) > 0:
                    key = list(depth_dict.values())[-1]
                    max_depth = max(depth_dict.keys())
                else:
                    key = "remainder"
                    max_depth = -1

            if name in split_strings or raw_rank == split_rank:
                key = name
                lines[key] = hierarchy.copy()
                depth_dict[depth] = name
                max_depth = depth
            elif add_hierarchy:
                for ancestor in hierarchy:
                    if ancestor not in lines[key]:
                        lines[key].append(ancestor)

            hierarchy.append(line)
            lines[key].append(line)

        for key in lines:
            save_file(key, lines[key])



# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        dest="report_file",
        required=True,
        help="Kraken or Bracken file of taxon relationships and quantities",
    )
    parser.add_argument(
        "--splits",
        dest="splits",
        required=False,
        nargs="*",
        default=[],
        help="List of taxon names to split the file by",
    )
    parser.add_argument(
            "--rank",
            dest="rank",
            required=False,
            help="The rank to split the file by",
        )

    args = parser.parse_args()

    rank_dict = {
        "kingdom": "K",
        "domain": "D",
        "phylum": "P",
        "class": "C",
        "order": "O",
        "family": "F",
        "genus": "G",
        "species": "S",
        "K": "K",
        "D": "D",
        "P": "P",
        "C": "C",
        "O": "O",
        "F": "F",
        "G": "G",
        "S": "S",
    }
    if args.rank:
        args.rank = rank_dict[args.rank]

    if not args.splits and not args.rank:
        args.splits=["Bacteria", "Viruses", "Metazoa"]



    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + "\n")

    parse_report_file(args.report_file, args.splits, args.rank)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
