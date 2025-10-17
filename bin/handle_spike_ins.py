#!/usr/bin/env python

import sys
import argparse
from datetime import datetime
import json
import os
from collections import defaultdict
import gzip


def parse_spike_in_refs(ref_file, ref_dir):
    ref_dict = {}
    with open(ref_file, "r") as f:
        ref_dict = json.load(f)
    if ref_dir:
        rel_path = os.path.abspath(ref_dir)
    else:
        ref_path = os.path.abspath(ref_file)
        rel_path = "/".join(os.path.abspath(ref_file).split("/")[:-1])
    for key in ref_dict:
        if ref_dict[key]["ref"] and not ref_dict[key]["ref"].startswith("/"):
            ref_dict[key]["ref"] = os.path.join(rel_path, ref_dict[key]["ref"])
        if ref_dict[key]["taxa"]:
            ref_dict[key]["taxa"] = [str(s) for s in ref_dict[key]["taxa"]]
    return ref_dict


def expand_spike_in_input(list_spike_ins, spike_in_dict):
    spike_taxids = []
    spike_refs = []
    for spike in list_spike_ins:
        if spike in spike_in_dict:
            sys.stdout.write(f"Named spike in {spike} found in reference dict\n")
            spike_taxids.extend(spike_in_dict[spike]["taxa"])
            if spike_in_dict[spike]["ref"]:
                spike_refs.append(spike_in_dict[spike]["ref"])
        elif spike.isnumeric():
            sys.stdout.write(f"Spike in {spike} given as additional taxid\n")
            spike_taxids.append(spike)
        elif spike.endswith("f*a") or spike.endswith("f*a.gz"):
            sys.stdout.write(f"Spike in {spike} given as additional reference file\n")
            spike_refs.append(spike)
        else:
            str_list_keys = ",".join(spike_in_dict.keys())
            sys.stdout.write(
                f"Spike in {spike} does not correspond to named spike_ins [{str_list_keys}] and is not a taxid or a reference fasta\n"
            )
            sys.exit(6)
    return spike_taxids, spike_refs


def parse_idxstats(idxstats_file):
    counts = defaultdict(int)
    total_unmapped = 0
    with open(idxstats_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            name, _, mapped, unmapped = parts
            mapped = int(mapped)
            unmapped = int(unmapped)
            if name == "*":
                total_unmapped = unmapped
                continue
            counts[name] = mapped
    counts["total"] = sum(counts.values()) + total_unmapped
    return counts


def read_reference_headers(spike_refs):
    ref_headers = {}
    for ref_path in spike_refs:
        headers = []
        with gzip.open(ref_path, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    headers.append(line[1:].strip().split()[0])
        ref_headers[ref_path] = headers
    return ref_headers


def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i != "":
            break
        depth += 1
    depth = int(depth / 2)
    return depth


def save_file(name, lines, header):
    if len(lines) == 0:
        return
    outfile = ".".join([name, "kreport_split.txt"])
    with open(outfile, "w") as f:
        if header:
            f.write(header)
        for line in lines:
            f.write(line)


def parse_report_file(report_file, spike_in, save_json):
    depth_dict = {}
    spike_in_entry = False
    max_depth = -1
    entries = {}
    header = None

    # parses a kraken or bracken file
    with open(report_file, "r") as f:
        for line in f:
            if line.startswith("% of Seqs"):
                header = line
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
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            depth = parse_depth(name)
            name = name.strip()
            rank = raw_rank[0]

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name,
                "taxid": ncbi,
                "is_spike_in": False,
            }

            while depth <= max_depth:
                spike_in_entry = False
                del depth_dict[max_depth]
                if len(depth_dict) > 0:
                    max_depth = max(depth_dict.keys())
                else:
                    max_depth = -1

            if ncbi in spike_in:
                depth_dict[depth] = name
                max_depth = depth
                spike_in_entry = True

            if spike_in_entry:
                entries[ncbi]["is_spike_in"] = True

    if save_json:
        with open(
            report_file.replace(".filtered", "").replace(".txt", ".json"), "w"
        ) as outfile:
            json.dump(entries, outfile, indent=4, sort_keys=False)

    spike_entries = defaultdict(lambda: defaultdict(str))
    for taxid in entries:
        if entries[taxid]["is_spike_in"]:
            spike_entries[taxid]["name"] = entries[taxid]["name"]
            spike_entries[taxid]["human_readable"] = entries[taxid]["name"]
            spike_entries[taxid]["taxid"] = entries[taxid]["taxid"]
            spike_entries[taxid]["classified_percentage"] = entries[taxid]["percentage"]
            spike_entries[taxid]["classified_count"] = entries[taxid][
                "count_descendants"
            ]
            spike_entries[taxid].setdefault("mapped_count", 0)
            spike_entries[taxid].setdefault("mapped_percentage", 0.0)
    return spike_entries


def combine_report_and_map_counts(
    list_spike_ins, spike_in_dict, report_entries, map_counts, ref_headers
):
    spike_summary = defaultdict(lambda: {})
    for spike in list_spike_ins:
        spike_dict = defaultdict(lambda: {})
        if spike in spike_in_dict:
            spike_name = spike
            if spike_in_dict[spike]["ref"]:
                for long_name in ref_headers[spike_in_dict[spike]["ref"]]:
                    name, taxid, taxon_name = long_name.split("|")
                    taxon_name = taxon_name.replace("_", " ")

                    entry = report_entries.get(taxid)
                    if entry:
                        spike_dict[name].update(entry)
                    else:
                        spike_dict[name]["taxid"] = taxid
                        spike_dict[name]["human_readable"] = taxon_name
                        spike_dict[name].setdefault("classified_count", 0)
                        spike_dict[name].setdefault("classified_percentage", 0.0)
                    spike_dict[name]["mapped_count"] = map_counts[long_name]
                    spike_dict[name]["mapped_percentage"] = (
                        float(map_counts[long_name]) / map_counts["total"] * 100
                    )
            else:
                for taxid in spike_in_dict[spike]["taxa"]:
                    entry = report_entries.get(taxid)
                    if entry:
                        spike_dict[taxid].update(entry)

        elif spike.isnumeric():
            entry = report_entries[spike]
            spike_name = entry["name"]
            spike_dict[spike].update(entry)
        elif spike.endswith("f*a") or spike.endswith("f*a.gz"):
            spike_name = spike.split("/")[-1].split(".")[0]
            for long_name in ref_headers[spike]:
                name, taxid, taxon_name = long_name.split("|")
                taxon_name = taxon_name.replace("_", " ")

                entry = report_entries.get(taxid)
                if entry:
                    spike_dict[name].update(entry)
                else:
                    spike_dict[name]["taxid"] = taxid
                    spike_dict[name]["human_readable"] = taxon_name
                    spike_dict[name].setdefault("classified_count", 0)
                    spike_dict[name].setdefault("classified_percentage", 0.0)
                spike_dict[name]["mapped_count"] = map_counts[long_name]
                spike_dict[name]["mapped_percentage"] = (
                    float(map_counts[long_name]) / map_counts["total"] * 100
                )
        spike_summary[spike_name].update(spike_dict)
    with open("spike_count_summary.json", "w") as outfile:
        json.dump(spike_summary, outfile, indent=4, sort_keys=False)
    return spike_summary


def check_spike_summary(spike_summary):
    check = {}
    for spike in spike_summary:
        found_any = False
        found_all = True

        found_classified_any = False
        found_mapped_any = False
        found_classified_all = True
        found_mapped_all = True

        for ref in spike_summary[spike]:
            found_classified = False
            found_mapped = False
            classified_count = spike_summary[spike][ref].get("classified_count")
            mapped_count = spike_summary[spike][ref].get("mapped_count")
            if classified_count and int(classified_count) > 0:
                found_classified = True
                found_classified_any = True
            else:
                found_classified_all = False

            if mapped_count and int(mapped_count) > 0:
                found_mapped = True
                found_mapped_any = True
            else:
                found_mapped_all = False

            if found_classified or found_mapped:
                found_any = True
            else:
                found_all = False

        if not found_any:
            found_all = False
        if not found_classified_any:
            found_classified_all = False
        if not found_mapped_any:
            found_mapped_all = False

        if found_all:
            check[spike] = "detected"
        elif found_any:
            check[spike] = "partial"
        else:
            check[spike] = "absent"

        if found_classified_all:
            print(f"Spike {spike} found all refs by classification")
        elif found_classified_any:
            print(f"Spike {spike} found some refs by classification")
        else:
            print(f"Spike {spike} absent by classification")

        if found_mapped_all:
            print(f"Spike {spike} found all refs by mapping")
        elif found_mapped_any:
            print(f"Spike {spike} found some refs by mapping")
        else:
            print(f"Spike {spike} absent by mapping")

    with open("spike_summary.json", "w") as outfile:
        json.dump(check, outfile, indent=4, sort_keys=False)
    return check


# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        dest="report_file",
        required=True,
        help="Kraken or Bracken file of taxon relationships and quantities",
    )
    parser.add_argument(
        "--idxstats",
        dest="idxstats_file",
        required=True,
        help="samtools idxstats output generated during spike removal",
    )
    parser.add_argument(
        "--spike_ins",
        dest="spike_ins",
        required=False,
        nargs="*",
        default=[],
        help="List of spike in names",
    )
    parser.add_argument(
        "--spike_in_dict",
        dest="spike_in_dict",
        required=True,
        help="JSON file specifying the known spike ins and their reference files. Assumes reference filepaths are either absolute or relative to the spike_in_ref_dir",
    )
    parser.add_argument(
        "--spike_in_ref_dir",
        dest="spike_in_ref_dir",
        help="Assumes reference filepaths are either absolute or relative to this filepath",
    )
    parser.add_argument(
        "--save_json",
        action="store_true",
        required=False,
        help="Save the kraken report in JSON format",
    )

    args = parser.parse_args()
    spike_ins = []
    for spike in args.spike_ins:
        spike_ins.extend(spike.split(","))

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + "\n")

    spike_in_dict = parse_spike_in_refs(args.spike_in_dict, args.spike_in_ref_dir)
    spike_taxids, spike_refs = expand_spike_in_input(spike_ins, spike_in_dict)

    spike_kraken_entries = parse_report_file(
        args.report_file, spike_taxids, args.save_json
    )

    map_counts = parse_idxstats(args.idxstats_file)
    ref_headers = read_reference_headers(spike_refs)

    if spike_ins:
        spike_summary = combine_report_and_map_counts(
            spike_ins, spike_in_dict, spike_kraken_entries, map_counts, ref_headers
        )
        check_spike_summary(spike_summary)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
