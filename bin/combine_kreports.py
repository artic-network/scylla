#!/usr/bin/env python

import sys
import gzip
from Bio import SeqIO
import argparse
from datetime import datetime

def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i!="":
            break
        depth += 1
    depth = int(depth/2)
    return depth

def get_kraken_hierarchy(kraken_file, entries={}, total=0):
    with open(kraken_file, "r") as f:
        hierarchy = []
        for line in f:
            if line.startswith("% of Seqs"):
                continue
            percentage, num_clade_root, num_direct, raw_rank, ncbi, name = line.strip().split("\t")
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            total += num_direct
            depth = parse_depth(name)
            name = name.strip()
            rank = raw_rank[0]
            hierarchy = hierarchy[:depth]
            hierarchy.append(ncbi)

            if ncbi in entries:
                entries[ncbi]["count"] += num_direct
                entries[ncbi]["count_descendants"] += num_clade_root
            else:
                entries[ncbi] = {"percentage": percentage, "count": num_direct, "count_descendants": num_clade_root,
            "raw_rank": raw_rank, "rank": rank, "depth": depth, "ncbi": ncbi, "name": name, "parents":[], "children":set()}

            if len(hierarchy) > 1:
                parent = hierarchy[-2]
                assert entries[parent]["depth"] < entries[ncbi]["depth"]
                entries[ncbi]["parents"] = hierarchy[:-1]
                entries[parent]["children"].add(ncbi)

    return entries,total

def combine_kraken_reports(kreports):
    entries = {}
    total = 0
    for kreport in kreports:
        entries, total = get_kraken_hierarchy(kreport, entries, total)
        print(kreport, len(entries), total)
    for ncbi in entries:
        entries[ncbi]["percentage"] = entries[ncbi]["count_descendants"] / float(total) *100
    return entries

def write_entry(out_handle, entry):
    offset = entry["depth"]
    out_handle.write("%f\t%i\t%i\t%s\t%s\t%s%s\n" %(entry["percentage"], entry["count"], entry["count_descendants"], entry["raw_rank"], entry["ncbi"], 2*offset*" ", entry["name"] ))

def choose_next(entries, taxa, ncbi="0"):
    if len(taxa) == 0:
        return None, taxa

    next = None
    if ncbi == "0":
        next = "1"
        taxa.remove(ncbi)
        taxa.remove(next)
    else:
        max_pcent = 0
        for child in entries[ncbi]["children"]:
            if child in taxa:
                if entries[child]["percentage"] >= max_pcent:
                    next = child
                    max_pcent = entries[child]["percentage"]
        if not next:
            parent = entries[ncbi]["parents"][-1]
            next, taxa = choose_next(entries, taxa, ncbi=parent)
        else:
            taxa.remove(next)

    return next, taxa

def write_kraken_report(outfile, entries):
    taxa = list(entries.keys())
    with open(outfile, "w") as out_report:
        out_report.write("% of Seqs\tClades\tTaxonomies\tRank\tTaxonomy ID\tScientific Name\n")
        current = "0"
        while len(taxa) > 0:
            write_entry(out_report, entries[current])
            current, taxa = choose_next(entries, taxa, ncbi=current)

def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='kraken_report_files', required=True, nargs='+',
        help='Kraken report files to parse')
    parser.add_argument('-o', "--outfile",dest='outfile', required=True,
        help='Output name')

    args=parser.parse_args()

    #Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')

    entries = combine_kraken_reports(args.kraken_report_files)
    write_kraken_report(args.outfile, entries)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
