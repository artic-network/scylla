#!/usr/bin/env python
"""Script to create an aggregated count from lineage data."""
### taken from epi2ome/wf-metagenomics under the conditions of their license https://github.com/epi2me-labs/wf-metagenomics/blob/master/LICENSE and modified for this use case

import argparse
import json
import sys

UNCLASSIFIED = 'Unclassified'
UNKNOWN = 'Unknown'

RANKS = [
    "superkingdom",
    "clade",
    "kingdom",
    "phylum",
    "subphylum"
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "serotype"
]


def update_or_create_unclassified(entries, unclassified_count):
    """Handle unclassified entries."""
    entries[UNCLASSIFIED] = {
        'rank': RANKS[0],
        'count': int(unclassified_count),
        'children': {
            UNKNOWN: {
                'rank': "species",
                'count': int(unclassified_count),
                'children': {}
            }
        }
    }
    return entries


def update_or_create_count(entry, entries, bracken_counts):
    """Increment lineage counts given entries."""
    tax_id, lineage, ranks = entry.rstrip().split('\t')
    lineage_split = lineage.split(';')
    ranks_split = ranks.split(';')
    count = int(bracken_counts[tax_id])

    previous = entries
    previous_rank = None
    for [name, rank] in zip(lineage_split, ranks_split):

        if rank not in RANKS:
            if previous_rank == "species":
                rank = "subspecies"
            else:
                continue

        current = previous.get(name)
        if not current:
            new_entry = {
                'rank': rank,
                'count': count,
                'children': {}
            }
            previous[name] = new_entry
            previous = new_entry['children']
            continue

        current['count'] += count
        previous = current['children']
        previous_rank = rank

    return entries


def yield_entries(entries, total, indent=0):
    """Get entries in printable form."""
    for i, j in entries.items():
        perc = "{:.2%}".format(j['count'] / total)
        yield (indent, i, j['count'], perc, j['rank'])
        for k in yield_entries(j['children'], total, indent + 1):
            yield k


def main(prefix, lineages, bracken, report):
    """Run lineage aggregation algorithm."""
    bracken_counts = {}
    entries = {}
    total = 0
    with open(bracken) as f:
        bracken = f.readlines()
    if len(bracken) > 0:
        for i in bracken:
            bracken_counts[i.split()[1]] = i.split()[0]
        with open(lineages) as f:
            infile = f.readlines()
        for line in infile:
            try:
                entries = update_or_create_count(line, entries, bracken_counts)
                total += 1
            except ValueError:
                sys.stderr.write(
                    """Lineage for tax id {} not found in taxonomy database"""
                    .format(str(line)))
    with open(report) as f:
        report_file = f.readlines()
        for line in report_file:
            if line.split()[4] == "0" and "unclassified" in line:
                unclassified_count = line.split()[2]
                entries = update_or_create_unclassified(
                    entries, unclassified_count)
                total += int(unclassified_count)
    output_report = open('{}.lineages.txt'.format(prefix), 'w')
    output_json = open('{}.lineages.json'.format(prefix), 'w')
    for entry in yield_entries(entries, total):
        [indent, name, count, perc, rank] = entry
        output_report.write(' '.join(
            ['-' * (indent + 1), name, str(count), perc, rank, '\n'])
        )
    output_json.write(json.dumps(entries))


def execute():
    """Parse command line arguments and run main."""
    parser = argparse.ArgumentParser(
        description="Aggregates lineage counts in a kraken2-like format",
    )
    parser.add_argument(
        '-i',
        help=(
            "Lineages .tsv (taxid, lineage)."
        ),
        dest="lineages",
    )

    parser.add_argument(
        '-b',
        help=(
            "Bracken Lineages .tsv (taxid, count)."
        ),
        dest="bracken",
    )

    parser.add_argument(
        '-u',
        help=(
            "full report to get unclassified count"
        ),
        dest="report",
    )

    parser.add_argument(
        '-p',
        help="Prefix to append to output file names.",
        dest="prefix",
        required=True,
    )

    args = parser.parse_args()

    main(
        lineages=args.lineages,
        prefix=args.prefix,
        bracken=args.bracken,
        report=args.report
    )


if __name__ == "__main__":
    execute()
