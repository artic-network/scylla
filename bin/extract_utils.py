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
import statistics as stats

from assignment import trim_read_id


def mean(l):
    """
    Takes a list of numbers and returns the mean.

    Args:
        l (list): A list of numbers.

    Returns:
        0: If l had zero length.
        float: The mean of items in the list.
    """
    if len(l) == 0:
        return 0
    return stats.fmean(l)


def median(l):
    """
    Takes a list of numbers and returns the median.

    Args:
        l (list): A list of numbers.

    Returns:
        0: If l had zero length.
        float: The median of items in the list.
    """
    if len(l) == 0:
        return 0
    return stats.median(l)


def check_read_files(reads):
    """
    Takes a read filename and checks if the file is zipped, a FASTA or FASTQ format.

    Args:
        reads (str): The reads filename.

    Returns:
        filetype (str): Either "fasta" or "fastq".
        zipped (bool): Whether the file is gzipped or not.
    """
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

def setup_outfiles(fastq_2, prefixes, filetype, get_handles=True):
    """
    Sets up the output filenames, optionally opening them and returning handles.

    Args:
        fastq_2 (bool): True or a string if fastq_2 exists, False or None if reads were not paired.
        prefixes (dict): A dictionary from the key (usually taxon_id) to the prefix of the output file for that key.
        filetype (str): "fasta" or "fastq" depending on input filetype.
        get_handles (bool): Whether to open file handles for the output files.

    Returns:
        filenames (dict): Dictionary with keys from the prefixes, to a list of output filenames.
        out_handles_1 (dict): Key (matching filenames) to out_handle (if opened)
        out_handles_2 (dict): Key (matching filenames) to out_handle (if opened and paired reads)
    """

    out_handles_1 = {}
    out_handles_2 = {}
    filenames = defaultdict(list)

    for key, prefix in prefixes.items():
        if fastq_2:
            filenames[key].append(f"{prefix}_1.{filetype}")
            filenames[key].append(f"{prefix}_2.{filetype}")
            if get_handles:
                if key in out_handles_1:
                    print("already open")
                else:
                    out_handles_1[key] = open(f"{prefix}_1.{filetype}", "w")
                    out_handles_2[key] = open(f"{prefix}_2.{filetype}", "w")

        else:
            filenames[key].append(f"{prefix}.{filetype}")
            if get_handles:
                if key in out_handles_1:
                    print("already open")
                else:
                    out_handles_1[key] = open(f"{prefix}.{filetype}", "w")
    return filenames, out_handles_1, out_handles_2

def close_outfiles(out_handles_1, out_handles_2):
    """
    Checks each handle in the dictionary and closes if it is open.

    Args:
        out_handles_1 (dict): Key to out_handle
        out_handles_1 (dict): Key to out_handle (if paired reads)
    """
    for handle in out_handles_1:
        out_handles_1[handle].close()
    for handle in out_handles_2:
        out_handles_2[handle].close()

def update_summary_with_record(taxon_id, record, out_counts, quals, lens):
    """
    Updates the summary dictionaries for counts, qualities and lengths of reads from the new record Dictionaries indexed by taxon.

    Args:
        taxon_id (str): Key for dictionaries
        record (SeqRecord): A read parsed by Bio.SeqIO
        out_counts (dict): Taxon ID to read count for that taxon.
        quals (dict): Taxon ID to list of quality scores (each averaged over the read) for that taxon.
        lens (dict): Taxon ID to list of read lengths for that taxon.
    """
    name, seq, qual = record
    out_counts[taxon_id] += 1
    quals[taxon_id].append(median([ord(x) - 33 for x in qual]))
    lens[taxon_id].append(len(seq))

def add_record(taxon_id, record, out_counts, quals, lens, filenames, file_index, out_handles={}, out_records={}, buffer_record_size=None):
    """
    Add the record to the required out file and update the summary statistics with it. If out_handles has open
    out_handles, reads will be written directly to the handle. If not it will be stored in memory in the out_records. If
    a buffer size is set and the out_records structure has at least this number of reads, the out_records for this taxon
    ID will be appended to the relevant file and cleared.

    Args:
        taxon_id (str): Key for dictionaries
        record (SeqRecord): A read parsed by Bio.SeqIO
        out_counts (dict): Taxon ID to read count for that taxon.
        quals (dict): Taxon ID to list of quality scores (each averaged over the read) for that taxon.
        lens (dict): Taxon ID to list of read lengths for that taxon.
        filenames (dict): Taxon ID to list of output filenames.
        file_index (int): 0 or 1, depending if processing forward or reverse read, index for the filenames list.
        out_handles (dict): {} if handles closed, Taxon ID to open handle if open.
        out_records (dict): {} if handles open, Taxon ID to list of records if not.
        buffer_record_size (int): Number of records per taxon ID before they should be output to reduce memory use.
    """

    name, seq, qual = record
    update_summary_with_record(taxon_id, record, out_counts, quals, lens)

    if out_handles != {} and taxon_id in out_handles:
        out_handles[taxon_id].write(f"@{name}\n{seq}\n+\n{qual}\n")
    else:
        if not buffer_record_size or len(out_records[taxon_id]) < buffer_record_size:
            out_records[taxon_id].append(record)
        else:
            with open(filenames[taxon_id][file_index], "a") as f:
                for existing_record in out_records[taxon_id]:
                    f.write(f"@{name}\n{seq}\n+\n{qual}\n")
            del out_records[taxon_id]
            out_records[taxon_id].append(record)

def save_records_to_file(out_records, filenames, file_index):
    """
    Save the records to file.

    Args:
        out_records (dict): Taxon ID to list of records
        filenames (dict): Taxon ID to list of output filenames
        file_index (int): 0 or 1, depending if processing forward or reverse read, index for the filenames list.
    """
    sys.stderr.write(f"Writing records for file {file_index+1}\n")
    for taxon, records in out_records.items():
        with open(filenames[taxon][file_index], "a") as f:
            for record in records:
                name, seq, qual = record
                f.write(f"@{name}\n{seq}\n+\n{qual}\n")
    del out_records

def file_iterator(read_file, read_map, subtaxa_map, inverse, out_counts, quals, lens, filenames, file_index, out_handles, low_memory=False):
    """
    Iterate through the read_file file and add the read to the appropriate file or handle.

    Args:
        read_file (str): Name of read file.
        read_map (dict): Dictionary from read ID to list of taxon ID associated with it.
        subtaxa_map (dict): Dictionary from taxon ID (output by read map) to list of taxon ID associated with
                            out files to be updated
        inverse (bool): If True, exclude all reads which match the taxon IDs in the read_map. If False, select them.
        out_counts (dict): Taxon ID to read count for that taxon.
        quals (dict): Taxon ID to list of quality scores (each averaged over the read) for that taxon.
        lens (dict): Taxon ID to list of read lengths for that taxon.
        filenames (dict): Taxon ID to list of output filenames
        file_index (int): 0 or 1, depending if processing forward or reverse read, index for the filenames list.
        out_handles (dict): {} if handles closed, Taxon ID to open handle if open.
        low_memory (bool): If True we expect to write directly to out handles, if False store records then write to file.
    Returns:
        int: Count of reads written to file.
    """
    reads_of_interest = set(read_map.keys())
    count = 0
    out_records = None
    if not low_memory:
        out_records = defaultdict(list)

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record
        trimmed_name = trim_read_id(name)
        if inverse:
            if trimmed_name in reads_of_interest:
                for taxon in read_map[trimmed_name]:
                    update_summary_with_record(taxon, record, out_counts, quals, lens)
                continue
            else:
                count += 1
                add_record("other", record, out_counts, quals, lens, filenames, file_index, out_handles=out_handles, out_records=out_records)

        elif not inverse:
            if trimmed_name not in reads_of_interest:
                continue
            for k2_taxon in read_map[trimmed_name]:
                for taxon in subtaxa_map[k2_taxon]:
                    count += 1
                    add_record(taxon, record, out_counts, quals, lens, filenames, file_index, out_handles=out_handles, out_records=out_records)

    if not low_memory:
        save_records_to_file(out_records, filenames, file_index)

    return count

def process_read_files(
    prefixes: str,
    filetype: str,
    read_map: dict,
    subtaxa_map: dict,
    read_file_1: Path,
    read_file_2: Path = None,
    inverse: bool = False,
    get_handles: bool = False
) -> tuple[dict, dict, dict, dict]:
    """
    Iterate through (paired/unpaired) read_files and save the relevant reads, collecting summary statistics for these reads.

    Args:
        prefixes (dict): A dictionary from the key (usually taxon_id) to the prefix of the output file for that key.
        filetype (str): "fasta" or "fastq" depending on input filetype.
        read_map (dict): Dictionary from read ID to list of taxon ID associated with it.
        subtaxa_map (dict): Dictionary from taxon ID (output by read map) to list of taxon ID associated with
                            out files to be updated.
        read_file_1 (Path): Name of read file.
        read_file_2 (Path): Name of read file pair (if it exists).
        inverse (bool): If True, exclude all reads which match the taxon IDs in the read_map. If False, select them.
        get_handles (bool): If True we expect to write directly to out handles, if False store records then write to file.

    Returns:
        out_counts (dict): Taxon ID to read count for that taxon.
        quals (dict): Taxon ID to list of quality scores (each averaged over the read) for that taxon.
        lens (dict): Taxon ID to list of read lengths for that taxon.
        filenames (dict): Taxon ID to list of output filenames.
    """
    out_counts = defaultdict(int)
    quals = defaultdict(list)
    lens = defaultdict(list)

    filenames, out_handles_1, out_handles_2 = setup_outfiles(fastq_2, prefixes, filetype, get_handles=get_handles)

    forward_count = file_iterator(read_file_1, read_map, subtaxa_map, inverse, out_counts, quals, lens, filenames, 0, out_handles_1, low_memory=get_handles)

    reverse_count = 0
    if fastq_2:
        reverse_count = file_iterator(read_file_2, read_map, subtaxa_map, inverse, out_counts, quals, lens, filenames, 1, out_handles_2, low_memory=get_handles)

        if forward_count != reverse_count and (forward_count == 0 or reverse_count == 0):
            sys.stderr.write(
                "ERROR: No reads found for one of the file pair: extracted %i an %i reads respectively" % (forward_count, reverse_count)
            )
            sys.exit(7)
    close_outfiles(out_handles_1, out_handles_2)
    return (out_counts, quals, lens, filenames)


def generate_summary(lists_to_extract, entries, prefix, out_counts, quals, lens, filenames, include_unclassified=False, short=False):
    """
    Generate a summary JSON file, including information about each taxon ID and corresponding read statistics.

    Args:
        lists_to_extract (list): A list of taxon ID to include in the summary.
        entries (dict): A dictionary containing information about the name/rank associated with each taxon ID, either
                        inferred from the Kraken report or NCBI taxonomy.
        prefix (str): Prefix for the summary file.
        out_counts (dict): Taxon ID to read count for that taxon.
        quals (dict): Taxon ID to list of quality scores (each averaged over the read) for that taxon.
        lens (dict): Taxon ID to list of read lengths for that taxon.
        filenames (dict): Taxon ID to list of output filenames.
        include_unclassified (bool): True if each output file includes unclassified reads as well as those associated
                                     with the taxon ID.
        short (bool): Output short form summary, excluding taxon information.
        filenames (dict): Taxon ID to list of output filenames.
    """
    sys.stderr.write("Write summary\n")
    summary = []
    for taxon_id in out_counts:
        if out_counts[taxon_id] == 0:
            sys.stderr.write(f"No reads extracted  for taxon_id {taxon_id}\n")
            continue

        if short:
            summary.append(
            {
                "filenames": filenames[taxon_id],
                "qc_metrics": {
                    "num_reads": out_counts[taxon_id],
                    "avg_qual": mean(quals[taxon_id]),
                    "mean_len": mean(lens[taxon_id]),
                },
                "includes_unclassified": include_unclassified
            }
            )
        else:
            summary.append(
            {
                "human_readable": entries[taxon_id].name,
                "taxon_id": taxon_id,
                "tax_level": entries[taxon_id].rank,
                "filenames": filenames[taxon_id],
                "qc_metrics": {
                    "num_reads": out_counts[taxon_id],
                    "avg_qual": mean(quals[taxon_id]),
                    "mean_len": mean(lens[taxon_id]),
                },
                "includes_unclassified": include_unclassified
            }
            )

    with open(f"{prefix}_summary.json", "w") as f:
        json.dump(summary, f)