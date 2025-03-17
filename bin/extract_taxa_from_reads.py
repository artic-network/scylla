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

from extract_utils import *
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
    """
    Loops through the kraken report, and for each taxon_id in the report, if it meets the thresholds, a key is added to
    lists_to_extract. The values in this dictionary are all taxon_ids which should be added to the file alongside the
    key taxon_id (e.g. parents or children)
    """
    lists_to_extract = defaultdict(set)
    for taxon in kraken_report.entries:
        entry = kraken_report.entries[taxon]
        if len(target_ranks) > 0 and entry.rank not in target_ranks:
            continue
        if min_count and entry.ucount < min_count:
            continue
        if min_count_descendants and entry.count < min_count_descendants:
            continue
        if min_percent and kraken_report.get_percentage(taxon, denominator=entry.domain) < min_percent:
            continue
        if len(names) > 0 and entry.name not in names and taxon not in names:
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


def setup_prefixes(lists_to_extract, prefix=None):
    outprefix = {}
    #if prefix:
    #    for taxid in lists_to_extract:
    #        outprefix[taxid] =  f"{prefix}_{taxid}"
    #else:
    for taxid in lists_to_extract:
        outprefix[taxid] =  f"{taxid}"
    return outprefix
  

def extract_taxa(
    kraken_report, lists_to_extract, kraken_assignment, reads1, reads2, prefix
):
    # open read files
    filetype, zipped = check_read_files(reads1)

    # This inverts the lists_to_extract to identify for an assigned taxon_id, which taxon_id files it should be added to.
    subtaxa_map = defaultdict(set)
    for taxon, subtaxons in lists_to_extract.items():
        for subtaxa in subtaxons:
            subtaxa_map[subtaxa].add(taxon)
        # sys.stderr.write(
        #    "INCLUDING PARENTS/CHILDREN, HAVE %i TAXA TO INCLUDE IN READ FILES for %s\n"
        #    % (len(lists_to_extract[taxon]), taxon)
        # )
    read_map = kraken_assignment.get_read_map(subtaxa_map)

    prefixes = setup_prefixes(lists_to_extract, prefix)
    out_counts, quals, lens, filenames, total_length = process_read_files(
        prefixes, filetype, read_map, subtaxa_map, reads1, reads2, inverse=False, get_handles=False
    )

    generate_summary(lists_to_extract, kraken_report.entries, prefix, out_counts, quals, lens, filenames, total_length, short=False)
    return out_counts

def check_out_counts(out_counts, kraken_report):
    for taxon_id in out_counts:
        if out_counts[taxon_id] != kraken_report.entries[taxon_id].count:
            sys.stderr.write(
                f"Did not correctly extract all reads for {taxon_id}, found {out_counts[taxon_id]} whilst the kraken report contains {kraken_report.entries[taxon_id].count}"
            )
            sys.exit(2)

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
        required=False,
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
        help="Minimum percentage of classified reads e.g 4",
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
    sys.stderr.write("Loading kraken report\n")
    kraken_report = KrakenReport(args.report_file)
    if args.max_human:
        kraken_report.check_host({"9606":args.max_human})

    # Initialize kraken assignment file
    sys.stderr.write("Loading kraken assignments\n")
    kraken_assignment = KrakenAssignments(args.kraken_assignment_file)

    sys.stderr.write("Identifying lists to extract\n")
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

    sys.stderr.write("Extracting reads from file\n")
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
