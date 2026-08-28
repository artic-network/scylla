#!/usr/bin/env python

import sys
import argparse
import json as json_module
from datetime import datetime

from extract_utils import *
from assignment import KrakenAssignments
from taxonomy import Taxonomy


def setup_prefixes(
    list_taxon_ids, entries, prefix, inverse=False, include_unclassified=False
):
    outprefix = {}
    if inverse:
        return {"other": prefix}

    for taxon_id in list_taxon_ids:
        taxon_name = entries[taxon_id].name.lower()
        if include_unclassified and taxon_name != "unclassified":
            taxon_name += "_and_unclassified"
        taxon_name = taxon_name.replace("viruses", "viral")
        outprefix[taxon_id] = taxon_name

    return outprefix


def build_fraction_spec(cfg, loaded_taxonomy, taxonomy_dir):
    """
    Build the (taxon_id_map, prefixes) pieces for one fraction config entry - the
    taxonomy-dependent setup that's identical whether this fraction is extracted alone or as part
    of a shared multi-fraction pass.

    Args:
        cfg (dict): A fraction config with keys 'prefix', 'taxid' (list), 'exclude' (bool) and
                    'include_unclassified' (bool).
        loaded_taxonomy (Taxonomy): A Taxonomy object with parent/child relationships already
                                     loaded. `entries` is populated (accumulated across calls) as
                                     needed by this function.
        taxonomy_dir (str): The unzipped directory downloaded from NCBI taxonomy.

    Returns:
        taxon_id_map (dict): Taxon ID (and descendants) to set of taxon_ids in cfg["taxid"].
        prefixes (dict): Key (usually taxon_id) to output filename prefix for this fraction.
    """
    taxon_id_map = loaded_taxonomy.get_taxon_id_map(
        cfg["taxid"], cfg["include_unclassified"]
    )
    # setup_prefixes only ever looks up entries for the originally-requested taxid (and not at
    # all when exclude is set, since it returns early), so loading taxonomy info for the whole
    # descendant closure in taxon_id_map - which can be the entire nodes.dmp file - is wasted work.
    if not cfg["exclude"]:
        loaded_taxonomy.load_entries(taxonomy_dir, cfg["taxid"])
    prefixes = setup_prefixes(
        cfg["taxid"],
        loaded_taxonomy.entries,
        cfg["prefix"],
        inverse=cfg["exclude"],
        include_unclassified=cfg["include_unclassified"],
    )
    return taxon_id_map, prefixes


def extract_fractions(fraction_configs, kraken_assignment, reads1, reads2, taxonomy_dir):
    """
    Extract one or more independent fractions in a single pass over the kraken assignment file
    and read file(s) - N separate fraction configs are not a partition of the read set (a read
    can legitimately belong to more than one fraction, e.g. a viral read is also "not human"), so
    this keeps each fraction's read_map/taxon_id_map/output files completely independent; only
    the (expensive, once-per-read) assignment file parse and FASTQ decompress+parse are shared.

    Args:
        fraction_configs (list): A list of fraction config dicts, each with keys 'prefix',
                                  'taxid' (list), 'exclude' (bool) and 'include_unclassified' (bool).
        kraken_assignment (KrakenAssignments): The loaded kraken assignment file.
        reads1 (str): FASTA/FASTQ file containing the raw reads.
        reads2 (str): 2nd FASTA/FASTQ file containing the raw reads (paired), or "" if unpaired.
        taxonomy_dir (str): The unzipped directory downloaded from NCBI taxonomy.

    Returns:
        dict: fraction prefix -> out_counts (taxon_id -> read count) for that fraction.
    """
    filetype, zipped = check_read_files(reads1)
    loaded_taxonomy = Taxonomy(taxonomy_dir)

    taxon_id_maps = []
    prefixes_list = []
    for cfg in fraction_configs:
        taxon_id_map, prefixes = build_fraction_spec(cfg, loaded_taxonomy, taxonomy_dir)
        taxon_id_maps.append(taxon_id_map)
        prefixes_list.append(prefixes)

    read_maps = kraken_assignment.get_read_maps(taxon_id_maps, paired=bool(reads2))

    fraction_specs = [
        {
            "prefixes": prefixes_list[i],
            "read_map": read_maps[i],
            "subtaxa_map": taxon_id_maps[i],
            "inverse": fraction_configs[i]["exclude"],
        }
        for i in range(len(fraction_configs))
    ]

    results, total_length = process_read_files_multi(
        fraction_specs, filetype, reads1, reads2
    )

    out_counts_by_fraction = {}
    for i, cfg in enumerate(fraction_configs):
        out_counts, read_stats, filenames = results[i]
        generate_summary(
            cfg["taxid"],
            loaded_taxonomy.entries,
            cfg["prefix"],
            out_counts,
            read_stats,
            filenames,
            total_length,
            include_unclassified=(cfg["include_unclassified"] != cfg["exclude"]),
            short=True,
            # total_length is the length of the whole input file - the same value for every
            # fraction here, so only the first fraction's call needs to write it out.
            write_total_length=(i == 0),
        )
        out_counts_by_fraction[cfg["prefix"]] = out_counts

    return out_counts_by_fraction


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
        required=False,
        help="Prefix for output files. Mutually exclusive with --fraction_config.",
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
        help="List of taxonomy ID[s] or names to exclude (space-delimited) from outputs",
    )
    parser.add_argument(
        "--fraction_config",
        dest="fraction_config",
        required=False,
        help=(
            "Path to a JSON file listing multiple fractions to extract in a single pass over the "
            "read file(s), instead of a single -p/--taxid/--exclude/--include_unclassified "
            "fraction. Each entry is an object with keys 'prefix', 'taxid' (list of taxonomy IDs "
            "or names), 'exclude' (bool) and 'include_unclassified' (bool). Mutually exclusive "
            "with -p/--prefix."
        ),
    )
    parser.set_defaults(append=False)

    args = parser.parse_args()

    if bool(args.prefix) == bool(args.fraction_config):
        parser.error(
            "Exactly one of -p/--prefix or --fraction_config must be provided."
        )

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    if args.fraction_config:
        with open(args.fraction_config) as f:
            fraction_configs = json_module.load(f)
    else:
        fraction_configs = [
            {
                "prefix": args.prefix,
                "taxid": args.taxid,
                "exclude": args.exclude,
                "include_unclassified": args.include_unclassified,
            }
        ]

    kraken_assignment = KrakenAssignments(args.kraken_assignment_file)
    out_counts_by_fraction = extract_fractions(
        fraction_configs, kraken_assignment, args.reads1, args.reads2, args.taxonomy
    )

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.stderr.write("READ COUNTS: \n")

    for prefix, out_counts in out_counts_by_fraction.items():
        for taxon in out_counts:
            sys.stderr.write("%s/%s: %i\n" % (prefix, taxon, out_counts[taxon]))

    sys.exit(0)


if __name__ == "__main__":
    main()
