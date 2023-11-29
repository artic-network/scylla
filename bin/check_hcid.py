#!/usr/bin/env python
import sys
import os
import argparse
import json
import csv
from datetime import datetime
from collections import defaultdict
import mappy as mp
import argparse
from Bio import SeqIO
import glob


def load_from_taxonomy(taxonomy_dir):
    taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
    parents = {}
    children = defaultdict(list)
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.split("\t|\t")
                tax_id, parent_tax_id = int(fields[0]), int(fields[1])
                parents[tax_id] = parent_tax_id
                children[parent_tax_id].append(tax_id)
    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy nodes.dmp file in %s" % taxonomy_dir
        )
        sys.exit(3)
    return parents, children

def parse_report_file(report_file, taxid_map):
    entries = {}
    counts = defaultdict(int)
    # parses a kraken or bracken file
    with open(report_file, "r") as f:
        for line in f:
            if line.startswith("% of Seqs"):
                continue
            try:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")
            except:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    ignore1,
                    ignore2,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")

            ncbi = int(ncbi)
            if ncbi not in taxid_map:
                continue
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            name = name.strip()
            rank = raw_rank[0]

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
            }
            counts[taxid_map[ncbi]]+=num_direct

    return entries, counts

def load_hcid_dict(hcid_file):
    hcid_dict = {}
    with open(hcid_file, "r") as f:
        hcid = json.load(f)
        for d in hcid:
            if d["taxon_id"]:
                hcid_dict[d["taxon_id"]] = d
                hcid_dict[d["taxon_id"]]["classified_count"] = 0
                hcid_dict[d["taxon_id"]]["classified_parent_count"] = 0
                hcid_dict[d["taxon_id"]]["classified_found"] = False
                hcid_dict[d["taxon_id"]]["classified_parent_found"] = False
    return hcid_dict

def check_report_for_hcid(hcid_dict, taxonomy_dir, kreport_file):
    parents, children = load_from_taxonomy(taxonomy_dir)

    taxid_map = {}
    for d in hcid_dict:
        if not d:
            continue
        taxids = set()
        lookup = [d]
        if "alt_taxon_ids" in hcid_dict:
            lookup.extend(hcid_dict["alt_taxon_ids"])
        while len(lookup) > 0:
            child = lookup.pop()
            if child in children:
                lookup.extend(children[child])
            taxid_map[child] = d
        if d in parents:
            taxid_map[parents[d]] = parents[d]

    entries, counts = parse_report_file(kreport_file, taxid_map)
    for taxid in hcid_dict:
        hcid_dict[taxid]["classified_count"] = counts[taxid]
        if taxid in parents and parents[taxid] in entries:
            print(taxid, parents[taxid], parents[taxid] in taxid_map, parents[taxid] in entries)
            hcid_dict[taxid]["classified_parent_count"]= entries[parents[taxid]]["count_descendants"]
        if counts[taxid] > hcid_dict[taxid]["min_count"]:
            hcid_dict[taxid]["classified_found"] = True
        if counts[parents[taxid]] > hcid_dict[taxid]["min_count"]:
            hcid_dict[taxid]["classified_parent_found"] = True


def map_to_refs(query, reference):
    counts = defaultdict(int)
    a = mp.Aligner(reference, best_n=1)  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")

    for name, seq, qual in mp.fastx_read(query): # read a fasta/q sequence
        for hit in a.map(seq): # traverse alignments
            counts[hit.ctg] += 1
            #print("{}\t{}\t{}\t{}\t{}".format(name, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            break
    return counts

def check_ref_coverage(hcid_dict, query, reference):
    counts = map_to_refs(query, reference)

    for taxon in hcid_dict:
        taxon_found = True
        hcid_dict[taxon]["mapped_count"] = 0
        hcid_dict[taxon]["mapped_required"] = 0
        hcid_dict[taxon]["mapped_required_details"] = []
        hcid_dict[taxon]["mapped_additional"] = 0
        hcid_dict[taxon]["mapped_found"] = True
        for ref in hcid_dict[taxon]["required_refs"]:
            if counts[ref] < hcid_dict[taxon]["min_count"]:
                hcid_dict[taxon]["mapped_found"] = False
            if counts[ref] > 0:
                hcid_dict[taxon]["mapped_required"] += 1
                hcid_dict[taxon]["mapped_count"] += counts[ref]
            hcid_dict[taxon]["mapped_required_details"].append("%s:%s" %(ref,str(counts[ref])))
        hcid_dict[taxon]["mapped_required"] = float(hcid_dict[taxon]["mapped_required"])/len(hcid_dict[taxon]["required_refs"])
        hcid_dict[taxon]["mapped_required_details"] = "|".join(hcid_dict[taxon]["mapped_required_details"])
        for ref in hcid_dict[taxon]["additional_refs"]:
            if counts[ref] > 0:
                hcid_dict[taxon]["mapped_additional"] += 1
                hcid_dict[taxon]["mapped_count"] += counts[ref]
        if len(hcid_dict[taxon]["additional_refs"])>0:
            hcid_dict[taxon]["mapped_additional"] = float(hcid_dict[taxon]["mapped_additional"])/len(hcid_dict[taxon]["additional_refs"])


def report_findings(hcid_dict, prefix):
    found = []
    for taxid in hcid_dict:
        if hcid_dict[taxid]["mapped_found"] and (hcid_dict[taxid]["classified_found"] or hcid_dict[taxid]["classified_parent_found"]) :
            with open("%s.warning" %taxid, "w") as f_warn:
                f_warn.write(f"WARNING: Found {hcid_dict[taxid]['classified_count']} classified reads ({hcid_dict[taxid]['mapped_count']} mapped reads) of {hcid_dict[taxid]['name']}\n")
            found.append(taxid)
    keys = ["name", "taxon_id", "min_count", "classified_count","mapped_count", "mapped_required", "mapped_required_details", "mapped_additional"]
    with open("%s.counts.csv" %prefix,'w') as f_counts:
        w = csv.DictWriter(f_counts,keys)
        w.writeheader()
        for taxid in hcid_dict:
            w.writerow({key: hcid_dict[taxid][key] for key in keys})
    return found

# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-k",
        dest="kreport_file",
        required=False,
        help="Kraken or Bracken file of taxon relationships and quantities",
    )
    parser.add_argument(
        "-r",
        dest="reads",
        required=True,
        help="FASTQ of reads",
    )
    parser.add_argument(
        "-t",
        dest="taxonomy",
        required=False,
        help="Taxonomy directory containing the nodes.dmp file. If not provided will infer from report file but this may lead to fewer reads extracted",
    )
    parser.add_argument(
        "-i",
        dest="hcid_file",
        required=False,
        default="resources/hcid.json",
        help="JSON file specifying HCID to look for a min counts to notify",
    )
    parser.add_argument(
        "-p",
        dest="prefix",
        required=False,
        default="hcid",
        help="Output prefix",
    )
    parser.add_argument(
        "-d",
        dest="ref_fasta",
        required=False,
        default="resources/hcid_refs.fq.gz",
        help="Reference FASTA for each HCID",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    hcid_dict = load_hcid_dict(args.hcid_file)
    check_report_for_hcid(hcid_dict, args.taxonomy, args.kreport_file)
    check_ref_coverage(hcid_dict, args.reads, args.ref_fasta)
    report_findings(hcid_dict, args.prefix)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
