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

from report import KrakenReport
from taxonomy import Taxonomy


def load_hcid_dict(hcid_file):
    hcid_dict = {}
    with open(hcid_file, "r") as f:
        hcid = json.load(f)
        for d in hcid:
            if d["taxon_id"]:
                d["taxon_id"] = str(d["taxon_id"])
                hcid_dict[d["taxon_id"]] = d
                hcid_dict[d["taxon_id"]]["classified_count"] = 0
                hcid_dict[d["taxon_id"]]["classified_parent_count"] = 0
                hcid_dict[d["taxon_id"]]["classified_found"] = False
                hcid_dict[d["taxon_id"]]["classified_parent_found"] = False
    return hcid_dict

def infer_taxid_map(hcid_dict, loaded_taxonomy):
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
            if child in loaded_taxonomy.children:
                lookup.extend(loaded_taxonomy.children[child])
            taxid_map[child] = d
        if d in loaded_taxonomy.parents:
            taxid_map[loaded_taxonomy.parents[d]] = loaded_taxonomy.parents[d]
    return taxid_map

def check_report_for_hcid(hcid_dict, taxonomy_dir, kreport_file):
    loaded_taxonomy = Taxonomy(taxonomy_dir)
    taxid_map = infer_taxid_map(hcid_dict, loaded_taxonomy)

    kraken_report = KrakenReport(kreport_file)
    for taxid in hcid_dict:
        hcid_dict[taxid]["classified_count"] = kraken_report.entries[taxid].count
        if taxid in loaded_taxonomy.parents and loaded_taxonomy.parents[taxid] in kraken_report.entries:
            print(
                taxid,
                loaded_taxonomy.parents[taxid],
                loaded_taxonomy.parents[taxid] in taxid_map,
                loaded_taxonomy.parents[taxid] in kraken_report.entries,
            )
            hcid_dict[taxid]["classified_parent_count"] = kraken_report.entries[loaded_taxonomy.parents[taxid]].count
        if hcid_dict[taxid]["classified_count"] > hcid_dict[taxid]["min_count"]:
            hcid_dict[taxid]["classified_found"] = True
        if hcid_dict[taxid]["classified_parent_count"] > hcid_dict[taxid]["min_count"]:
            hcid_dict[taxid]["classified_parent_found"] = True


def map_to_refs(query, reference, preset):
    counts = defaultdict(int)
    ranges = defaultdict(list)
    a = mp.Aligner(reference, best_n=1, preset=preset)  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")

    read_count = 0
    for name, seq, qual in mp.fastx_read(query):  # read a fasta/q sequence
        read_count += 1
        for hit in a.map(seq):  # traverse alignments
            counts[hit.ctg] += 1
            ranges[hit.ctg].append([hit.r_st, hit.r_en])
            # print("{}\t{}\t{}\t{}\t{}".format(name, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            break
        #if read_count % 1000000 == 0:
        #    break
    return counts, ranges


def check_pileup(ref, ref_ranges, reference_file, min_coverage=0):
    if len(ref_ranges) == 0:
        return 0
    for name, seq, qual in mp.fastx_read(reference_file):
        if name == ref:
            coverages = [0] * len(seq)
            for r in ref_ranges:
                for i in range(r[0], r[1]):
                    coverages[i] += 1
            zeros = [i for i in coverages if i <= min_coverage]
            return float(len(seq) - len(zeros)) / len(seq)
    return 0


def check_ref_coverage(hcid_dict, query, reference, preset):
    counts, ranges = map_to_refs(query, reference, preset)

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
            ref_covg = check_pileup(ref, ranges[ref], reference)
            if ref_covg < hcid_dict[taxon]["min_coverage"]:
                hcid_dict[taxon]["mapped_found"] = False
            hcid_dict[taxon]["mapped_required_details"].append(
                "%s:%i:%f" % (ref, counts[ref], ref_covg)
            )
        hcid_dict[taxon]["mapped_required"] = float(
            hcid_dict[taxon]["mapped_required"]
        ) / len(hcid_dict[taxon]["required_refs"])
        hcid_dict[taxon]["mapped_required_details"] = "|".join(
            hcid_dict[taxon]["mapped_required_details"]
        )
        for ref in hcid_dict[taxon]["additional_refs"]:
            if counts[ref] > 0:
                hcid_dict[taxon]["mapped_additional"] += 1
                hcid_dict[taxon]["mapped_count"] += counts[ref]
        if len(hcid_dict[taxon]["additional_refs"]) > 0:
            hcid_dict[taxon]["mapped_additional"] = float(
                hcid_dict[taxon]["mapped_additional"]
            ) / len(hcid_dict[taxon]["additional_refs"])


def report_findings(hcid_dict, prefix):
    found = []
    for taxid in hcid_dict:
        if hcid_dict[taxid]["mapped_found"] and (
            hcid_dict[taxid]["classified_found"]
            or hcid_dict[taxid]["classified_parent_found"]
        ):
            with open("%s.warning.json" % taxid, "w") as f_warn:
                msg1 = f"WARNING: Found {hcid_dict[taxid]['classified_count']} classified reads ({hcid_dict[taxid]['mapped_count']} mapped reads) of {hcid_dict[taxid]['name']} and {hcid_dict[taxid]['classified_parent_count']} classified reads for the parent taxon.\n"
                msg2 = f"Mapping details for required references (ref_accession:mapped_read_count:fraction_ref_covered) {hcid_dict[taxid]['mapped_required_details']}.\n"
                warning =   {
                                "msg":msg1+msg2,
                                "taxid":taxid,
                                "classified_count":hcid_dict[taxid]['classified_count'],
                                "mapped_count":hcid_dict[taxid]['mapped_count'],
                                "mapped_details":f"ref_accession:mapped_read_count:fraction_ref_covered|{hcid_dict[taxid]['mapped_required_details']}"
                            }
                json.dump(warning, f_warn, indent=4, sort_keys=False)
            found.append(taxid)
    keys = [
        "name",
        "taxon_id",
        "min_count",
        "classified_count",
        "classified_parent_count",
        "mapped_count",
        "mapped_required",
        "mapped_required_details",
        "mapped_additional",
    ]
    with open("%s.counts.csv" % prefix, "w") as f_counts:
        w = csv.DictWriter(f_counts, keys)
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
        default="resources/hcid_refs.fa.gz",
        help="Reference FASTA for each HCID",
    )
    parser.add_argument(
            "--illumina",
            action="store_true",
            required=False,
            help="Use the short read minimap preset",
        )

    args = parser.parse_args()
    preset = None
    if args.illumina:
        preset = "sr"

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    sys.stderr.write("Load HCID information\n")
    hcid_dict = load_hcid_dict(args.hcid_file)
    sys.stderr.write("Check kraken report for counts\n")
    check_report_for_hcid(hcid_dict, args.taxonomy, args.kreport_file)
    sys.stderr.write("Check mapped coverage\n")
    check_ref_coverage(hcid_dict, args.reads, args.ref_fasta, preset)
    sys.stderr.write("Report findings\n")
    report_findings(hcid_dict, args.prefix)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
