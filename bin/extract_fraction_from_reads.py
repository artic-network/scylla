#!/usr/bin/env python

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

def setup_prefixes(list_taxon_ids, entries, prefix, inverse=False, include_unclassified=False):
    outprefix = {}
    if inverse:
        return {"other": prefix}

    for taxon_id in list_taxon_ids:
        taxon_name = entries[taxon_id].name.lower()
        if include_unclassified and taxon_name != "unclassified":
            taxon_name += "_and_unclassified"
        taxon_name = taxon_name.replace("viruses", "viral")
        outprefix[taxon_id] = taxon_name

    return outprefix
  

def extract_reads(
    read_map, taxon_id_map, entries, reads1, reads2, prefix, taxon_ids, exclude, include_unclassified
):
    # check read files
    filetype, zipped = check_read_files(reads1)

    prefixes = setup_prefixes(taxon_ids, entries, prefix, inverse=exclude, include_unclassified=include_unclassified)
    out_counts, quals, lens, filenames = process_read_files(prefixes, filetype, read_map, taxon_id_map, reads1, reads2, inverse=exclude, get_handles=True)

    generate_summary(taxon_ids, entries, prefix, out_counts, quals, lens, filenames, include_unclassified=(include_unclassified != exclude), short=True)

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
        "--include_unclassified",
        dest="include_unclassified",
        required=False,
        action="store_true",
        default=False,
        help="Include unclassified in output files",
    )
    parser.add_argument(
        "--exclude",
        dest="exclude",
        action="store_true",
        default=False,
        help="List of taxonomy ID[s] or names to exclude (space-delimited) from outputs"
    )
    parser.set_defaults(append=False)

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    loaded_taxonomy = Taxonomy(args.taxonomy)
    taxon_id_map = loaded_taxonomy.get_taxon_id_map(args.taxid, args.include_unclassified)
    loaded_taxonomy.load_entries(args.taxonomy, taxon_id_map.keys())

    # Initialize kraken assignment file
    kraken_assignment = KrakenAssignments(args.kraken_assignment_file)
    read_map = kraken_assignment.parse_kraken_assignment_file(taxon_id_map)

    out_counts = extract_reads(
        read_map,
        taxon_id_map,
        loaded_taxonomy.entries,
        args.reads1,
        args.reads2,
        args.prefix,
        args.taxid,
        args.exclude,
        args.include_unclassified
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
