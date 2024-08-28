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

from extract_utils import mean,median,check_read_files
from report import KrakenReport
from assignment import KrakenAssignments
from taxonomy import Taxonomy

def get_taxon_id_lists(
    kraken_report,
    loaded_taxonomy,
    names=[],
    target_ranks=[],
    min_count=None,
    min_count_descendants=None,
    min_percent=None,
    top_n=None,
    include_parents=False,
    include_children=False
):
    lists_to_extract = defaultdict(set)
    for taxon in kraken_report:
        entry = kraken_report[taxon]
        if len(target_ranks) > 0 and entry["rank"] not in target_ranks:
            continue
        if min_count and entry["ucount"] < min_count:
            continue
        if min_count_descendants and entry["count"] < min_count_descendants:
            continue
        if min_percent and kraken_report.get_percentage(taxon) < min_percent:
            continue
        if len(names) > 0 and entry["name"] not in names and taxon not in names:
            continue

        lists_to_extract[taxon].add(taxon)
        if include_parents:
            lookup = [taxon]
            while len(lookup) > 0:
                parent = lookup.pop()
                if parent in loaded_taxonomy.parents and parent != "1":
                    lookup.append(loaded_taxonomy.parents[lookup])
                if parent in kraken_report.entries and parent != "1":
                    lookup.append(kraken_report.entries[parent].parent)
                if parent != "1":
                    lists_to_extract[taxon].add(parent)

        if include_children:
            lookup = [taxon]
            while len(lookup) > 0:
                child = lookup.pop()
                lists_to_extract[taxon].add(child)
                if child in loaded_taxonomy.children:
                    lookup.extend(loaded_taxonomy.children[child])
                if child in kraken_report.entries:
                    lookup.extend(kraken_report.entries[child].children)
    sys.stderr.write("SELECTED %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    if top_n and len(lists_to_extract) > top_n:
        X = list(lists_to_extract.keys())
        Y = [kraken_report.get_percentage(x) for x in X]
        ordered = [x for _, x in sorted(zip(Y, X))]
        to_delete = ordered[top_n:]
        for taxon in to_delete:
            del lists_to_extract[taxon]
        sys.stderr.write("REDUCED TO %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    return lists_to_extract

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

    out_records_1 = defaultdict(list)
    out_records_2 = defaultdict(list)

    sys.stderr.write("Reading in\n")
    for record in pyfastx.Fastq(fastq_1, build_index=False):
        name, seq, qual = record
        trimmed_name = trim_read_id(name)
        if trimmed_name not in reads_of_interest:
            continue

        for k2_taxon in read_map[trimmed_name]:
            for taxon in subtaxa_map[k2_taxon]:
                out_counts[taxon] += 1
                quals[taxon].append(median([ord(x) - 33 for x in qual]))
                lens[taxon].append(len(seq))

                out_records_1[taxon].append(record)

    sys.stderr.write("Writing records\n")
    for taxon, records in out_records_1.items():
        if fastq_2:
            with open(f"{taxon}_1.{filetype}", "w") as f:
                for record in records:
                    name, seq, qual = record
                    f.write(f"@{name}\n{seq}\n+\n{qual}\n")
        else:
            with open(f"{taxon}.{filetype}", "w") as f:
                for record in records:
                    name, seq, qual = record
                    f.write(f"@{name}\n{seq}\n+\n{qual}\n")
    del out_records_1

    if fastq_2:
        sys.stderr.write("Reading in second file of pair\n")
        for record in pyfastx.Fastq(fastq_2, build_index=False):
            name, seq, qual = record
            trimmed_name = trim_read_id(name)
            if trimmed_name not in reads_of_interest:
                continue

            for k2_taxon in read_map[trimmed_name]:
                for taxon in subtaxa_map[k2_taxon]:
                    out_counts[taxon] += 1
                    quals[taxon].append(median([ord(x) - 33 for x in qual]))
                    lens[taxon].append(len(seq))

                    out_records_2[taxon].append(record)

        sys.stderr.write("Writing records for second file in pair\n")
        for taxon, records in out_records_2.items():
            with open(f"{taxon}_2.{filetype}", "w") as f:
                for record in records:
                    name, seq, qual = record
                    f.write(f"@{name}\n{seq}\n+\n{qual}\n")

    return (out_counts, quals, lens)

def extract_taxa(
    kraken_report, lists_to_extract, kraken_assignment, reads1, reads2, prefix
):
    # open read files
    filetype, zipped = check_read_files(reads1)

    subtaxa_map = defaultdict(set)

    for taxon, subtaxons in lists_to_extract.items():
        for subtaxa in subtaxons:
            subtaxa_map[subtaxa].add(taxon)

        # sys.stderr.write(
        #    "INCLUDING PARENTS/CHILDREN, HAVE %i TAXA TO INCLUDE IN READ FILES for %s\n"
        #    % (len(lists_to_extract[taxon]), taxon)
        # )
    read_map = kraken_assignment.parse_kraken_assignment_file(subtaxa_map)

    sys.stderr.write("Iterating through read file\n")
    out_counts, quals, lens = fastq_iterator(
        prefix, filetype, read_map, subtaxa_map, reads1, reads2
    )

    sys.stderr.write("Write summary\n")
    summary = []
    for taxon in lists_to_extract:
        if out_counts[taxon] == 0:
            sys.stderr.write("No reads extracted  for taxid %s\n" %taxon)
            continue
        if reads2:
            summary.append(
                {
                    "human_readable": kraken_report[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": kraken_report[taxon]["rank"],
                    "filenames": [
                        "%s_1.%s" % (taxon, filetype),
                        "%s_2.%s" % (taxon, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    }
                }
            )
        else:
            summary.append(
                {
                    "human_readable": kraken_report[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": kraken_report[taxon]["rank"],
                    "filenames": [
                        "%s.%s" % (taxon, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    }
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
        required=True,
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

    loaded_taxonomy = None
    if args.taxonomy:
        sys.stderr.write("Loading taxonomy\n")
        loaded_taxonomy = Taxonomy(args.taxonomy)

    # Load kraken report entries
    sys.stderr.write("Loading kreport\n")
    kraken_report = KrakenReport(args.report_file)
    kraken_report.check_host({9606:args.max_human})

    # Initialize kraken assignment file
    kraken_assignment = KrakenAssignment(args.kraken_assignment_file)

    sys.stderr.write("Checking for lists to extract\n")
    lists_to_extract = get_taxon_id_lists(
        kraken_report,
        loaded_taxonomy,
        names=args.taxid,
        target_ranks=target_ranks,
        min_count=args.min_count,
        min_count_descendants=args.min_count_descendants,
        min_percent=args.min_percent,
        top_n=args.top_n,
        include_parents=args.include_parents,
        include_children=args.include_children
    )

    sys.stderr.write("Performing extractions\n")
    out_counts = extract_taxa(
        kraken_report,
        lists_to_extract,
        kraken_assignment,
        args.reads1,
        args.reads2,
        args.prefix
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
