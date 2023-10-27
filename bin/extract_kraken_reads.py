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
import pyfastx
import argparse
import json
from datetime import datetime
from collections import defaultdict
from pathlib import Path


def mean(l):
    if len(l) == 0:
        return 0
    return sum(l) / len(l)


def median(l):
    if len(l) % 2 == 0:
        i = (len(l)) / 2
    else:
        i = (len(l) + 1) / 2
    i = int(i)
    l = sorted(l)
    return l[i]


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
        sys.exit(3)
    return parents, children


def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i != "":
            break
        depth += 1
    depth = int(depth / 2)
    return depth


def infer_hierarchy(report_file):
    parents = {}
    children = defaultdict(list)
    hierarchy = []
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
            hierarchy = hierarchy[: depth - 1]
            hierarchy.append(ncbi)

            if len(hierarchy) > 1:
                parent = hierarchy[-2]
                parents[ncbi] = parent
                children[parent].append(ncbi)
    return parents, children


def load_report_file(report_file, max_human=None):
    entries = {}
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
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            name = name.strip()
            rank = raw_rank[0]

            if name in ["Homo sapiens", "unclassified", "root"]:
                if max_human and name == "Homo sapiens" and num_direct > max_human:
                    sys.stderr.write(
                        "ERROR: found %i human reads, max allowed is %i\n"
                        % (num_direct, max_human)
                    )
                    sys.exit(2)
                continue

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
            }

    sys.stderr.write("FOUND %i TAXA IN BRACKEN REPORT\n" % len(entries))
    return entries


def get_taxon_id_lists(
    report_entries,
    parents,
    children,
    names=[],
    target_ranks=[],
    min_count=None,
    min_count_descendants=None,
    min_percent=None,
    top_n=None,
    include_parents=False,
    include_children=False,
):
    lists_to_extract = {}
    for taxon in report_entries:
        entry = report_entries[taxon]
        if len(target_ranks) > 0 and entry["rank"] not in target_ranks:
            continue
        if min_count and entry["count"] < min_count:
            continue
        if min_count_descendants and entry["count_descendants"] < min_count_descendants:
            continue
        if min_percent and entry["percentage"] < min_percent:
            continue
        if len(names) > 0 and entry["name"] not in names and taxon not in names:
            continue

        lists_to_extract[taxon] = [taxon]
        if include_parents:
            lookup = taxon
            while lookup in parents and lookup != "1":
                lookup = parents[lookup]
                if lookup != "1":
                    lists_to_extract[taxon].append(lookup)

        if include_children:
            lookup = [taxon]
            while len(lookup) > 0:
                child = lookup.pop()
                if child != taxon:
                    lists_to_extract[taxon].append(child)
                lookup.extend(children[child])
    sys.stderr.write("SELECTED %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    if top_n and len(lists_to_extract) > top_n:
        X = list(lists_to_extract.keys())
        Y = [report_entries[x]["percentage"] for x in X]
        ordered = [x for _, x in sorted(zip(Y, X))]
        to_delete = ordered[top_n:]
        for taxon in to_delete:
            del lists_to_extract[taxon]
        sys.stderr.write("REDUCED TO %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    return lists_to_extract


def check_read_files(reads):
    if reads[-3:] == ".gz":
        read_file = gzip.open(reads, "rt")
        zipped = True
    else:
        read_file = open(reads, "rt")
        zipped = False
    first = read_file.readline()
    if len(first) == 0:
        sys.stderr.write("ERROR: sequence file's first line is blank\n")
        sys.exit(1)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(1)
    return filetype, zipped


def parse_kraken_assignment_line(line):
    line_vals = line.strip().split("\t")
    if len(line_vals) < 5:
        return -1, ""
    if "taxid" in line_vals[2]:
        temp = line_vals[2].split("taxid ")[-1]
        tax_id = temp[:-1]
    else:
        tax_id = line_vals[2]

    read_id = trim_read_id(line_vals[1])

    if tax_id == "A":
        tax_id = 81077
    else:
        tax_id = tax_id
    return tax_id, read_id


def trim_read_id(read_id):
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


def fastq_iterator(
    prefix: str,
    filetype: str,
    read_map: dict,
    subtaxa_map: dict,
    fastq_1: Path,
    fastq_2: Path = None,
) -> tuple[dict, dict, dict]:
    """Func to iterate over fastq files and extract reads of interest

    Args:
        prefix (str): Outfile prefix
        filetype (str): output filetype (only affects naming)
        read_map (dict): dict of read_id: taxon_list (from kraken report)
        subtaxa_map (dict): dict of subtaxa: output taxon (from lists_to_extract)
        fastq_1 (Path): Path to fastq _1 file
        fastq_2 (Path, optional): Path to fastq _2 file if input is paired data. Defaults to None.

    Returns:
        tuple[dict, dict, dict]: number of reads written by taxa, quality scores by taxa, sequence length by taxa
    """
    reads_of_interest = set(read_map.keys())

    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)

    out_records = {}

    fq_1 = pyfastx.Fastq(fastq_1)
    if fastq_2:
        fq_2 = pyfastx.Fastq(fastq_2)

        for record_1, record_2 in zip(fq_1, fq_2):
            trimmed_name = trim_read_id(record_1.name)
            if trimmed_name not in reads_of_interest:
                continue

            if trimmed_name != trim_read_id(record_2.name):
                sys.stderr.write("ERROR: read names do not match\n")
                sys.exit(10)

            for k2_taxon in read_map[trimmed_name]:
                for taxon in subtaxa_map[k2_taxon]:
                    out_counts[taxon] += 2
                    quals[taxon].extend(
                        [median(record_1.quali), median(record_2.quali)]
                    )
                    lens[taxon].extend([len(record_1.seq), len(record_2.seq)])

                    try:
                        out_records[taxon][1].append(record_1)
                        out_records[taxon][2].append(record_2)

                    except KeyError:
                        out_records[taxon] = {1: [], 2: []}

                        out_records[taxon][1].append(record_1)
                        out_records[taxon][2].append(record_2)

        for taxon, records in out_records.items():
            with open(f"{prefix}.{taxon}_1.{filetype}", "w") as f:
                for record in records[1]:
                    f.write(record.raw)
            with open(f"{prefix}.{taxon}_2.{filetype}", "w") as f:
                for record in records[2]:
                    f.write(record.raw)

    else:
        for record in fq_1:
            trimmed_name = trim_read_id(record.name)
            if trimmed_name not in reads_of_interest:
                continue

            for k2_taxon in read_map[trimmed_name]:
                for taxon in subtaxa_map[k2_taxon]:
                    out_counts[taxon] += 1
                    quals[taxon].append(median(record.quali))
                    lens[taxon].append(len(record.seq))

                    try:
                        out_records[taxon].append(record)

                    except KeyError:
                        out_records[taxon] = []

                        out_records[taxon].append(record)

        for taxon, records in out_records.items():
            with open(f"{prefix}.{taxon}.{filetype}", "w") as f:
                for record in records:
                    f.write(record.raw)

    return (out_counts, quals, lens)


def extract_taxa(
    report_entries, lists_to_extract, kraken_assignment_file, reads1, reads2, prefix
):
    # open read files
    filetype, zipped = check_read_files(reads1)

    read_map = defaultdict(set)

    subtaxa_map = defaultdict(set)

    for taxon, subtaxons in lists_to_extract.items():
        for subtaxa in subtaxons:
            subtaxa_map[subtaxa].add(taxon)

        # sys.stderr.write(
        #    "INCLUDING PARENTS/CHILDREN, HAVE %i TAXA TO INCLUDE IN READ FILES for %s\n"
        #    % (len(lists_to_extract[taxon]), taxon)
        # )

    with open(kraken_assignment_file, "r") as kfile:
        for line in kfile:
            tax_id, read_id = parse_kraken_assignment_line(line)
            if tax_id in subtaxa_map.keys():
                read_map[read_id].add(tax_id)

    out_counts, quals, lens = fastq_iterator(
        prefix, filetype, read_map, subtaxa_map, reads1, reads2
    )

    summary = []
    for taxon in lists_to_extract:
        if reads2:
            summary.append(
                {
                    "human_readable": report_entries[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": report_entries[taxon]["rank"],
                    "filenames": [
                        "%s.%s_1.%s" % (prefix, taxon, filetype),
                        "%s.%s_2.%s" % (prefix, taxon, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                }
            )
        else:
            summary.append(
                {
                    "human_readable": report_entries[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": report_entries[taxon]["rank"],
                    "filenames": [
                        "%s.%s.%s" % (prefix, taxon, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                }
            )
    with open("%s_summary.json" % prefix, "w") as f:
        json.dump(summary, f)
    return out_counts


# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-k",
        dest="kraken_assignment_file",
        required=True,
        help="Kraken assignment file to parse",
    )
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

    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        required=True,
        default="taxid",
        help="Prefix for output files",
    )

    parser.add_argument(
        "--taxid",
        dest="taxid",
        required=False,
        nargs="*",
        default=[],
        help="List of taxonomy ID[s] or names to extract (space-delimited) - each to their own file",
    )
    parser.add_argument(
        "--rank", dest="rank", required=False, nargs="*", help="Rank(s) to extract"
    )
    parser.add_argument(
        "--max_human",
        dest="max_human",
        required=False,
        type=int,
        help="Maximum human reads to allow",
    )
    parser.add_argument(
        "--min_count",
        dest="min_count",
        required=False,
        type=int,
        help="Minimum direct read count",
    )
    parser.add_argument(
        "--min_count_descendants",
        dest="min_count_descendants",
        required=False,
        type=int,
        help="Minimum read count at taxon level or descendants",
    )
    parser.add_argument(
        "--min_percent",
        dest="min_percent",
        required=False,
        type=float,
        help="Minimum percentage of reads e.g 4",
    )
    parser.add_argument(
        "--n",
        dest="top_n",
        required=False,
        type=int,
        help="Maximum number of taxa to extract (top n)",
    )
    parser.add_argument(
        "--include_parents",
        dest="include_parents",
        required=False,
        action="store_true",
        default=False,
        help="Include reads classified at parent levels of the specified taxids",
    )
    parser.add_argument(
        "--include_children",
        dest="include_children",
        required=False,
        action="store_true",
        default=False,
        help="Include reads classified more specifically than the specified taxids",
    )
    parser.set_defaults(append=False)

    args = parser.parse_args()

    if not args.report_file:
        sys.stderr.write(
            "ERROR: require at least one report file from bracken or kraken\n"
        )
        sys.exit(1)

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

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
        target_ranks = [rank_dict[r] for r in args.rank]
    else:
        target_ranks = []

    parent, children = None, None
    if args.taxonomy:
        parent, children = load_from_taxonomy(args.taxonomy)
    else:
        parent, children = infer_hierarchy(args.report_file)

    # get taxids to extract
    report_entries = load_report_file(args.report_file, args.max_human)
    lists_to_extract = get_taxon_id_lists(
        report_entries,
        parent,
        children,
        names=args.taxid,
        target_ranks=target_ranks,
        min_count=args.min_count,
        min_count_descendants=args.min_count_descendants,
        min_percent=args.min_percent,
        top_n=args.top_n,
        include_parents=args.include_parents,
        include_children=args.include_children,
    )

    out_counts = extract_taxa(
        report_entries,
        lists_to_extract,
        args.kraken_assignment_file,
        args.reads1,
        args.reads2,
        args.prefix,
    )

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.stderr.write("READ COUNTS: \n")

    for taxon in out_counts:
        sys.stderr.write("%s: %i\n" % (taxon, out_counts[taxon]))

    sys.exit(0)


if __name__ == "__main__":
    main()
