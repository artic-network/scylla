#!/usr/bin/env python

import sys
import argparse
from datetime import datetime
import json


def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i != "":
            break
        depth += 1
    depth = int(depth / 2)
    return depth


def save_file(name, lines, header):
    if len(lines) == 0:
        return
    outfile = ".".join([name, "kreport_split.txt"])
    with open(outfile, "w") as f:
        if header:
            f.write(header)
        for line in lines:
            f.write(line)


def parse_report_file(report_file, split_strings, split_rank, ignore, save_json):
    depth_dict = {}
    lines = {"remainder": []}
    key = "remainder"
    ignore_entry = False
    hierarchy = []
    max_depth = -1
    entries = {}
    header = None

    # parses a kraken or bracken file
    with open(report_file, "r") as f:
        for line in f:
            if line.startswith("% of Seqs"):
                header = line
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
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            depth = parse_depth(name)
            name = name.strip()
            rank = raw_rank[0]

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
                "taxid": ncbi,
                "is_spike_in": False
            }

            hierarchy = hierarchy[:depth]
            add_hierarchy = False

            while depth <= max_depth:
                add_hierarchy = True
                ignore_entry = False
                del depth_dict[max_depth]
                if len(depth_dict) > 0:
                    key = list(depth_dict.values())[-1]
                    max_depth = max(depth_dict.keys())
                else:
                    key = "remainder"
                    max_depth = -1

            if ncbi in ignore:
                key = name
                depth_dict[depth] = name
                max_depth = depth
                ignore_entry = True
                add_hierarchy = False
            elif name in split_strings or raw_rank == split_rank:
                key = name
                lines[key] = hierarchy.copy()
                depth_dict[depth] = name
                max_depth = depth
                ignore_entry = False
            elif add_hierarchy:
                for ancestor in hierarchy:
                    if ancestor not in lines[key]:
                        lines[key].append(ancestor)

            hierarchy.append(line)
            if ignore_entry:
                entries[ncbi]["is_spike_in"] = True
            else:
                lines[key].append(line)

    for key in lines:
        save_file(key, lines[key], header)

    if save_json:
        with open(report_file.replace(".filtered", "").replace(".txt", ".json"), "w") as outfile:
            json.dump(entries, outfile, indent=4, sort_keys=False)

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
            "--ignore",
            dest="ignore",
            required=False,
            nargs="*",
            default=[],
            help="List of taxon ids to be ignored (along with their descendants)",
        )
    parser.add_argument(
        "--rank",
        dest="rank",
        required=False,
        help="The rank to split the file by",
    )
    parser.add_argument(
        "--save_json",
        action="store_true",
        required=False,
        help="Save the kraken report in JSON format",
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

    if args.ignore and len(args.ignore) > 0:
        sys.stdout.write("Filtering taxids: [%s]\n" %",".join(args.ignore))

    if not args.splits and not args.rank and not args.ignore:
        args.splits = ["Bacteria", "Viruses", "Metazoa"]

    if args.splits and len(args.splits) > 0:
        sys.stdout.write("Splitting on: [%s]\n" %",".join(args.splits))

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + "\n")

    parse_report_file(args.report_file, args.splits, args.rank, args.ignore, args.save_json)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
