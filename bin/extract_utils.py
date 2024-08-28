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
