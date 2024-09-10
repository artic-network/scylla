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

def setup_outfiles(fastq_2, list_prefix, filetype, open=True):
    out_handles_1 = {}
    out_handles_2 = {}
    filenames = defaultdict(list)

    for prefix in list_prefix:
        if fastq_2:
            filenames[prefix].append(f"{prefix}_1.{filetype}")
            filenames[prefix].append(f"{prefix}_2.{filetype}")
            if open:
                out_handles_1[prefix] = open(f"{prefix}_1.{filetype}", "w")
                out_handles_2[prefix] = open(f"{prefix}_2.{filetype}", "w")
        else:
            filenames[prefix].append(f"{prefix}.{filetype}")
            if open:
                out_handles_1[prefix] = open(f"{prefix}.{filetype}", "w")
    return filenames, out_handles_1, out_handles_2

def close_outfiles(out_handles_1, out_handles_2):
    for handle in out_handles_1:
        out_handles_1[handle].close()
    for handle in out_handles_2:
        out_handles_2[handle].close()

def update_summary_with_record(taxon, record, out_counts, quals, lens):
    name, seq, qual = record
    out_counts[taxon] += 1
    quals[taxon].append(median([ord(x) - 33 for x in qual]))
    lens[taxon].append(len(seq))

def add_record(taxon, record, out_counts, quals, lens, filenames, file_index, out_handles=None, out_records=None, buffer_record_size=None):
    if out_handles != {} and taxon in out_handles:
        update_summary_with_record(taxon, record, out_counts, quals, lens)
        out_handles[taxon].write(f"@{name}\n{seq}\n+\n{qual}\n")
    elif out_records:
        if not buffer_record_size or len(out_records[taxon]) < buffer_record_size:
            out_records[taxon].append(record)
        else:
            with open(filenames[taxon][file_index], "a") as f:
                for existing_record in out_records[taxon]:
                    name, seq, qual = existing_record
                    update_summary_with_record(taxon, existing_record, out_counts, quals, lens)
                    f.write(f"@{name}\n{seq}\n+\n{qual}\n")
            del out_records[taxon]
            out_records[taxon].append(record)

def save_records_to_file(out_records, filenames, file_index):
    sys.stderr.write("Writing records\n")
    for taxon, records in out_records.items():
        with open(filenames[taxon][file_index], "w") as f:
            for record in records:
                name, seq, qual = record
                f.write(f"@{name}\n{seq}\n+\n{qual}\n")
    del out_records

def file_iterator(fastq, read_map, subtaxa_map, inverse, file_index, out_handles, out_counts, quals, lens, filenames, low_memory=False):
    reads_of_interest = set(read_map.keys())
    out_records = None
    if not low_memory:
        out_records = defaultdict(list)

    sys.stderr.write("Reading in\n")
    for record in pyfastx.Fastq(fastq, build_index=False):
        name, seq, qual = record
        trimmed_name = trim_read_id(name)
        if inverse and trimmed_name in reads_of_interest:
            continue
        elif not inverse and trimmed_name not in reads_of_interest:
            continue

        for k2_taxon in read_map[trimmed_name]:
            for taxon in subtaxa_map[k2_taxon]:
                add_record(taxon, record, out_counts, quals, lens, filenames, out_handles=out_handles, out_records=out_records)

    if not low_memory:
        save_records_to_file(out_records, filenames, file_index)

    return out_records

def fastq_iterator(
    inverse: bool,
    filetype: str,
    read_map: dict,
    subtaxa_map: dict,
    fastq_1: Path,
    fastq_2: Path = None,
) -> tuple[dict, dict, dict, dict]:
    """Func to iterate over fastq files and extract reads of interest

    Args:
        filetype (str): output filetype (only affects naming)
        read_map (dict): dict of read_id: taxon_list (from kraken report)
        subtaxa_map (dict): dict of subtaxa: output taxon (from lists_to_extract)
        lists_to_extract (iter): list of taxon to extract a file for
        inverse (bool): should the reads kept be inverted (keep only those not in list)
        fastq_1 (Path): Path to fastq _1 file
        fastq_2 (Path, optional): Path to fastq _2 file if input is paired data. Defaults to None.

    Returns:
        tuple[dict, dict, dict]: number of reads written by taxa, quality scores by taxa, sequence length by taxa
    """
    reads_of_interest = set(read_map.keys())

    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)
    filenames, out_handles_1, out_handles_2 = setup_outfiles(fastq_1, fastq_2, list_prefix, filetype, open=False)

    out_records_1 = file_iterator(fastq_1, read_map, subtaxa_map, inverse, 0, out_handles_1, out_counts, quals, lens, filenames, low_memory=False):
    save_records_to_file(out_records_1, filenames, 0)

    if fastq_2:
        out_records_2 = file_iterator(fastq_2, read_map, subtaxa_map, inverse, 1, out_handles_2, out_counts, quals, lens, filenames, low_memory=False):
        save_records_to_file(out_records_2, filenames, 1)

    return (out_counts, quals, lens, filenames)


def generate_summary(lists_to_extract, entries, prefix, out_counts, quals, lens, filenames, includes_unclassified=False):
    sys.stderr.write("Write summary\n")
    summary = []
    for taxon in lists_to_extract:
        if out_counts[taxon] == 0:
            sys.stderr.write("No reads extracted  for taxid %s\n" %taxon)
            continue

        summary.append(
            {
                "human_readable": entries[taxon]["name"],
                "taxon": taxon,
                "tax_level": entries[taxon]["rank"],
                "filenames": filenames[taxon],
                "qc_metrics": {
                    "num_reads": out_counts[taxon],
                    "avg_qual": mean(quals[taxon]),
                    "mean_len": mean(lens[taxon]),
                }
                "includes_unclassified": includes_unclassified
            }
        )

    with open("%s_summary.json" % prefix, "w") as f:
        json.dump(summary, f)