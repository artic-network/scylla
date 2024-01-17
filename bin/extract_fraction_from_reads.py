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


def load_from_taxonomy(taxonomy_dir, taxids, include_unclassified):
    sys.stderr.write("Loading hierarchy\n")
    entries = {"0":{"name": "unclassified",
                    "taxon": "0",
                    "rank": None
                    }
               }


    taxid_map = defaultdict(set)
    for key in taxids:
        taxid_map[key].add(key)
    if include_unclassified:
        taxid_map["0"].update(taxids)

    children = defaultdict(set)
    parent = defaultdict(str)
    if len(taxids) > 0:
        try:
            taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
            with open(taxonomy, "r") as f:
                for line in f:
                    fields = line.split("\t|\t")
                    taxid, parent_taxid, rank = fields[0], fields[1], fields[2]
                    children[parent_taxid].add(taxid)
                    parent[taxid] = parent_taxid
                    if taxid in taxids:
                        entries[taxid] = {"name": None,
                                          "taxon": taxid,
                                          "rank": rank
                                          }
        except:
            sys.stderr.write(
                "ERROR: Could not find taxonomy nodes.dmp file in %s" % taxonomy_dir
            )
            sys.exit(4)

        try:
            taxonomy = os.path.join(taxonomy_dir, "names.dmp")
            with open(taxonomy, "r") as f:
                for line in f:
                    fields = line.split("\t|\t")
                    taxid, name, name_type = fields[0], fields[1], fields[3]
                    if taxid in taxids and ("scientific name" in name_type or entries[taxid]["name"] == None):
                        entries[taxid]["name"] = name

        except:
            sys.stderr.write(
                "ERROR: Could not find taxonomy names.dmp file in %s" % taxonomy_dir
            )
            sys.exit(4)

        check = list(taxids)
        while len(check) > 0:
            current = check.pop()
            check.extend(children[current])
            for child in children[current]:
                taxid_map[child].update(taxid_map[current])

    return taxid_map, entries, parent


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
        sys.exit(5)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(5)
    return filetype, zipped


def parse_kraken_assignment_line(line):
    line_vals = line.strip().split("\t")
    if len(line_vals) < 5:
        return -1, ""
    if "taxid" in line_vals[2]:
        temp = line_vals[2].split("taxid ")[-1]
        taxid = temp[:-1]
    else:
        taxid = line_vals[2]

    read_id = trim_read_id(line_vals[1])

    if taxid == "A":
        taxid = 81077
    else:
        taxid = taxid
    return taxid, read_id

def parse_kraken_assignment_file(kraken_assignment_file, taxid_map, parent):
    sys.stderr.write("Loading read assignments\n")
    read_map = defaultdict(set)
    with open(kraken_assignment_file, "r") as kfile:
        for line in kfile:
            taxid, read_id = parse_kraken_assignment_line(line)
            if taxid in taxid_map:
                read_map[read_id].update(taxid_map[taxid])
            else:
                # handle case where taxid has changed
                current = taxid
                while current in parent and current not in taxid_map and current != "1":
                    current = parent[current]
                    if current in taxid_map:
                        print(f"Add {taxid} to {current} list")
                        read_map[read_id].update(taxid_map[current])
            if read_id == "08c2e096-4393-f2b0-b327-155f13f52ecc":
                print(line)
                print(read_map[read_id])
    return read_map

def trim_read_id(read_id):
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


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
    names = defaultdict(list)

    sys.stderr.write("Open file handles\n")
    out_handles_1 = {}
    out_handles_2 = {}

    print(entries)
    for taxon, entry in entries.items():
        taxon_name = entry["name"].lower()
        if include_unclassified and taxon_name != "unclassified":
            taxon_name += "_and_unclassified"
        if fastq_2:
            out_handles_1[taxon] = open(f"{taxon_name}_1.{filetype}", "w")
            out_handles_2[taxon] = open(f"{taxon_name}_2.{filetype}", "w")
            names[taxon].append(f"{taxon_name}_1.{filetype}")
            names[taxon].append(f"{taxon_name}_2.{filetype}")
        else:
            out_handles_1[taxon] = open(f"{taxon_name}.{filetype}", "w")
            names[taxon].append(f"{taxon_name}.{filetype}")

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

    for taxon in out_handles_1:
        out_handles_1[taxon].close()
    for taxon in out_handles_2:
        out_handles_2[taxon].close()

    return (out_counts, quals, lens, names)


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

    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)
    names = defaultdict(list)

    sys.stderr.write("Open file handles\n")
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

    sys.stderr.write("Write summary\n")
    summary = []
    for taxon in names:
        if reads2:
            summary.append(
                {
                    "filenames": names[taxon],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                    "includes_unclassified": (include_unclassified != exclude) # this is xor
                }
            )
        else:
            summary.append(
                {
                    "filenames": names[taxon],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                    "includes_unclassified": (include_unclassified != exclude) # this is xor
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

    taxid_map, entries, parent = load_from_taxonomy(args.taxonomy, args.taxid, args.include_unclassified)
    read_map = parse_kraken_assignment_file(args.kraken_assignment_file, taxid_map, parent)

    out_counts = extract_reads(
        read_map,
        entries,
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
