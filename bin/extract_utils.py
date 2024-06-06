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
    if len(l) == 0:
        return 0
    if len(l) % 2 == 0:
        i = (len(l)) / 2
    else:
        i = (len(l) + 1) / 2
    i = int(i) - 1
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


def parse_kraken_assignment_file(kraken_assignment_file, taxid_map, parent=None):
    sys.stderr.write("Loading read assignments\n")
    read_map = defaultdict(set)
    with open(kraken_assignment_file, "r") as kfile:
        for line in kfile:
            taxid, read_id = parse_kraken_assignment_line(line)
            if taxid in taxid_map:
                read_map[read_id].update(taxid_map[taxid])
            elif parent:
                # handle case where taxid has changed
                current = taxid
                while current in parent and current not in taxid_map and current != "1":
                    current = parent[current]
                    if current in taxid_map:
                        print(f"Add {taxid} to {current} list")
                        read_map[read_id].update(taxid_map[current])
    return read_map


def trim_read_id(read_id):
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id
