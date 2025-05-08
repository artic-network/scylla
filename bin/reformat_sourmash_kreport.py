#!/usr/bin/env python

import sys
import os
import argparse
from datetime import datetime
from collections import defaultdict


def translate_names(taxonomy_dir, taxon_names):
    taxon_ids = defaultdict(str)
    name_dict = defaultdict(str)
    taxon_ids["unclassified"] = "0"
    name_dict["0"] = "unclassified"

    taxonomy = os.path.join(taxonomy_dir, "names.dmp")
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.strip().split("\t|")
                taxon_id, name, name_type = (
                    fields[0].strip(),
                    fields[1].strip(),
                    fields[3].strip(),
                )
                if name not in taxon_names and taxon_id not in taxon_names:
                    continue
                if name_type not in ["scientific name", "equivalent name"]:
                    continue
                taxon_ids[name] = taxon_id
                if name_type == "scientific name" or taxon_id not in name_dict:
                    name_dict[taxon_id] = name

    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy names.dmp file in %s" % taxonomy_dir
        )
        sys.exit(4)
    return taxon_ids, name_dict


def load_from_taxonomy(taxonomy_dir):
    taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
    parents = {}
    children = defaultdict(list)
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.split("\t|\t")
                tax_id, parent_tax_id = fields[0], fields[1]
                parents[tax_id] = parent_tax_id
                children[parent_tax_id].append(tax_id)
    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy nodes.dmp file in %s" % taxonomy_dir
        )
        sys.exit(4)
    return parents, children


def load_report_file(report_file):
    entries = defaultdict(lambda: defaultdict(str))
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
            name = name.strip()
            if name == "":
                continue
            rank = raw_rank[0]

            entries[name] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
                "taxid": None,
                "depth": 0,
            }

    sys.stdout.write("FOUND %i TAXA IN KREPORT FILE\n" % len(entries))
    return entries


def update_entries(entries, taxid_dict):
    for name in entries:
        if name in taxid_dict:
            entries[name]["taxid"] = taxid_dict[name]


def infer_order(entries, parents, names_dict):
    start = [entries[name]["taxid"] for name in entries if entries[name]["taxid"]]
    order = []
    while len(start) > 0:
        id = start.pop()
        if id == "0":
            order.insert(0, id)
            continue
        greedy = [id]
        while id in parents and parents[id] != "1":
            id = parents[id]
            if id in start:
                greedy.append(id)
                start.remove(id)
            if id in order:
                index = order.index(id) + 1
                depth = entries[names_dict[id]]["depth"]
                for i, taxid in enumerate(reversed(greedy)):
                    name = names_dict[taxid]
                    entries[name]["depth"] = i + depth + 1
                while len(greedy) > 0:
                    order.insert(index, greedy.pop(0))
                break

        order.extend(reversed(greedy))
        for i, taxid in enumerate(reversed(greedy)):
            name = names_dict[taxid]
            entries[name]["depth"] = i + 1

    return order


def write_kreport(outfile, entries, name_dict, order):
    with open(outfile, "w") as kreport:
        for taxid in order:
            name = name_dict[taxid]
            out_name = "  " * entries[name]["depth"] + name
            fields = [
                entries[name]["percentage"],
                entries[name]["count_descendants"],
                entries[name]["count"],
                entries[name]["raw_rank"],
                taxid,
                out_name,
            ]
            try:
                kreport.write("%s\n" % "\t".join(fields))
            except:
                print(fields)


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
        "-t",
        dest="taxonomy",
        required=True,
        help="Taxonomy directory containing the nodes.dmp file. If not provided will infer from report file but this may lead to fewer reads extracted",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="Output kreport file",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + "\n")

    entries = load_report_file(args.report_file)
    taxid_dict, name_dict = translate_names(args.taxonomy, entries.keys())
    update_entries(entries, taxid_dict)

    parents, children = load_from_taxonomy(args.taxonomy)

    order = infer_order(entries, parents, name_dict)
    write_kreport(args.outfile, entries, name_dict, order)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
