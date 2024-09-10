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

from extract_utils import mean,median,check_read_files
from report import KrakenReport
from assignment import KrakenAssignments
from taxonomy import Taxonomy


def fastq_iterator(
    prefix: str,
    filetype: str,
    include_unclassified: bool,
    entries: dict,
    read_map: dict,
    fastq_1: Path,
    fastq_2: Path = None,
) -> tuple[dict, dict, dict]:
    """Func to iterate over fastq files and extract reads of interest

    Args:
        prefix (str): Outfile prefix
        filetype (str): output filetype (only affects naming)
        include_unclassified (bool): true if includes unclassified reads
        read_map (dict): dict of read_id: taxon_list (from kraken report)
        subtaxa_map (dict): dict of subtaxa: output taxon (from load_from_taxonomy)
        exclude (bool): if true, inverts the logic for including reads
        fastq_1 (Path): Path to fastq _1 file
        fastq_2 (Path, optional): Path to fastq _2 file if input is paired data. Defaults to None.

    Returns:
        tuple[dict, dict, dict]: number of reads written by taxa, quality scores by taxa, sequence length by taxa
    """
    reads_of_interest = set(read_map.keys())

    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)

    print(entries)
    for taxon, entry in entries.items():
        taxon_name = entry["name"].lower()
        if include_unclassified and taxon_name != "unclassified":
            taxon_name += "_and_unclassified"
        taxon_name = taxon_name.replace("viruses", "viral")
        subtaxa_map[taxon_name] = [taxon]

    filenames, out_handles_1, out_handles_2 = setup_outfiles(fastq_2, subtaxa_map, filetype, open=True)


    sys.stderr.write("Iterating through read file\n")
    count = 0
    for record in pyfastx.Fastq(fastq_1, build_index=False):
        count += 1
        if count % 1000000 == 0:
            print(count)
        name, seq, qual = record
        trimmed_name = trim_read_id(name)
        if trimmed_name not in reads_of_interest:
            continue

        for taxon in read_map[trimmed_name]:
            out_handles_1[taxon].write(f"@{name}\n{seq}\n+\n{qual}\n")
            out_counts[taxon] += 1
            quals[taxon].append(median([ord(x) - 33 for x in qual]))
            lens[taxon].append(len(seq))

    if fastq_2:
        sys.stderr.write("Iterating second read file of pair\n")
        for record in pyfastx.Fastq(fastq_2, build_index=False):
            name, seq, qual = record
            trimmed_name = trim_read_id(name)
            if trimmed_name not in reads_of_interest:
                continue

            for taxon in read_map[trimmed_name]:
                out_handles_2[taxon].write(f"@{name}\n{seq}\n+\n{qual}\n")
                out_counts[taxon] += 1
                quals[taxon].append(median([ord(x) - 33 for x in qual]))
                lens[taxon].append(len(seq))

    close_outfiles(out_handles_1, out_handles_2)

    return (out_counts, quals, lens, names, filenames)



def fastq_iterator_inverse(
    prefix: str,
    filetype: str,
    taxids: list,
    read_map: dict,
    fastq_1: Path,
    fastq_2: Path = None,
) -> tuple[dict, dict, dict]:
    """Func to iterate over fastq files and extract reads of interest

    Args:
        prefix (str): Outfile prefix
        filetype (str): output filetype (only affects naming)
        read_map (dict): dict of read_id: taxon_list (from kraken report)
        subtaxa_map (dict): dict of subtaxa: output taxon (from load_from_taxonomy)
        exclude (bool): if true, inverts the logic for including reads
        fastq_1 (Path): Path to fastq _1 file
        fastq_2 (Path, optional): Path to fastq _2 file if input is paired data. Defaults to None.

    Returns:
        tuple[dict, dict, dict]: number of reads written by taxa, quality scores by taxa, sequence length by taxa
    """
    reads_of_interest = set(read_map.keys())
    print(reads_of_interest)

    filenames, out_handles_1, out_handles_2 =

    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)
    names = defaultdict(list)

    sys.stderr.write("Setup file handles\n")
    out_handles_1 = {}
    out_handles_2 = {}

    if fastq_2:
        out_handles_1["all"] = open(f"{prefix}_1.{filetype}", "w")
        out_handles_2["all"] = open(f"{prefix}_2.{filetype}", "w")
        names["all"].append(f"{prefix}_1.{filetype}")
        names["all"].append(f"{prefix}_2.{filetype}")
    else:
        out_handles_1["all"] = open(f"{prefix}.{filetype}", "w")
        names["all"].append(f"{prefix}.{filetype}")

    sys.stderr.write("Iterating through read file\n")
    count = 0
    for record in pyfastx.Fastq(fastq_1, build_index=False):
        count += 1
        if count % 1000000 == 0:
            print(count)
        name, seq, qual = record
        trimmed_name = trim_read_id(name)
        if trimmed_name  in reads_of_interest:
            continue

        out_handles_1["all"].write(f"@{name}\n{seq}\n+\n{qual}\n")
        out_counts["all"] += 1
        quals["all"].append(median([ord(x) - 33 for x in qual]))
        lens["all"].append(len(seq))

    if fastq_2:
        sys.stderr.write("Iterating second read file of pair\n")
        for record in pyfastx.Fastq(fastq_2, build_index=False):
            name, seq, qual = record
            trimmed_name = trim_read_id(name)
            if trimmed_name in reads_of_interest:
                continue

            out_handles_2["all"].write(f"@{name}\n{seq}\n+\n{qual}\n")
            out_counts["all"] += 1
            quals["all"].append(median([ord(x) - 33 for x in qual]))
            lens["all"].append(len(seq))

    for taxon in out_handles_1:
        out_handles_1[taxon].close()
    for taxon in out_handles_2:
        out_handles_2[taxon].close()

    return (out_counts, quals, lens, names)

def extract_reads(
    read_map, entries, reads1, reads2, prefix, taxids, exclude, include_unclassified
):
    # open read files
    filetype, zipped = check_read_files(reads1)

    if exclude:
        out_counts, quals, lens, names = fastq_iterator_inverse(
            prefix, filetype, taxids, read_map, reads1, reads2
        )
    else:
        out_counts, quals, lens, names = fastq_iterator(
            prefix, filetype, include_unclassified, entries, read_map, reads1, reads2
        )

    generate_summary(lists_to_extract, entries, out_counts, quals, lens, filenames, includes_unclassified=(include_unclassified != exclude))

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

    # Initialize kraken assignment file
    kraken_assignment = KrakenAssignment(args.kraken_assignment_file)
    read_map = kraken_assignment.parse_kraken_assignment_file(taxon_id_map)

    out_counts = extract_reads(
        read_map,
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
