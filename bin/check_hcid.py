#!/usr/bin/env python
import sys
import os
import argparse
import json
from datetime import datetime
from collections import defaultdict
from pathlib import Path


def load_from_taxonomy(taxonomy_dir):
    taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
    parents = {}
    children = defaultdict(list)
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.split("\t|\t")
                tax_id, parent_tax_id = int(fields[0]), int(fields[1])
                parents[tax_id] = parent_tax_id
                children[parent_tax_id].append(tax_id)
    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy nodes.dmp file in %s" % taxonomy_dir
        )
        sys.exit(3)
    return parents, children

def parse_report_file(report_file, taxid_map):
    entries = {}
    counts = defaultdict(int)
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

            ncbi = int(ncbi)
            if ncbi not in taxid_map:
                continue
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            name = name.strip()
            rank = raw_rank[0]

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
            }
            counts[taxid_map[ncbi]]+=num_direct

    return entries, counts

def check_report_for_hcid(taxonomy_dir, hcid_file, kreport_file):
    parents, children = load_from_taxonomy(taxonomy_dir)

    hcid_dict = {}
    with open(hcid_file, "r") as f:
        hcid = json.load(f)
        for d in hcid:
            if d["taxon_id"]:
                hcid_dict[d["taxon_id"]] = d

    taxid_map = {}
    for d in hcid_dict:
        if not d:
            continue
        taxids = set()
        lookup = [d]
        while len(lookup) > 0:
            child = lookup.pop()
            if child in children:
                lookup.extend(children[child])
            taxid_map[child] = d

    entries, counts = parse_report_file(kreport_file, taxid_map)
    for taxid in counts:
        hcid_dict[taxid]["found"] = counts[taxid]
        if counts[taxid] > hcid_dict[taxid]["count"]:
            sys.stdout.write(f"WARNING: Found {counts[taxid]} reads of {hcid_dict[taxid]['name']}\n")
    return hcid_dict

# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r",
        dest="report_file",
        required=False,
        help="Kraken or Bracken file of taxon relationships and quantities",
    )
    parser.add_argument(
        "-t",
        dest="taxonomy",
        required=False,
        help="Taxonomy directory containing the nodes.dmp file. If not provided will infer from report file but this may lead to fewer reads extracted",
    )
    parser.add_argument(
        "-i",
        dest="hcid_file",
        required=False,
        help="JSON file specifying HCID to look for a min counts to notify",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    check_report_for_hcid(args.taxonomy, args.hcid_file, args.report_file)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
