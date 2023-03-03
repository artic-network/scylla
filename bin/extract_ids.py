#!/usr/bin/env python

import argparse

def get_top_n_kraken_hits(braken_file,target_rank=None,min_num_reads=10,top_n=None,min_percent=None):
    with open(braken_file, "r") as f:
        percentage_dict = {}
        for line in f:
            percentage, num_clade_root, num_direct, rank, ncbi, name = line.strip().split("\t")
            name = name.strip()
            percentage = float(percentage)

            if name in ['Homo sapiens', 'unclassified', 'root']:
                continue
            elif name == 'unclassified':
                percentage_dict[ncbi] = percentage

            if target_rank and not rank.startswith(target_rank):
                continue
            if min_percent and percentage < min_percent:
                continue

            if int(num_direct) > min_num_reads:
                percentage_dict[ncbi] = percentage

    percentage_dict = dict(sorted(percentage_dict.items(), key=lambda item: item[1], reverse=True))
    if top_n:
        return list(percentage_dict)[:top_n]
    else:
        return list(percentage_dict)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--bracken_report", help="bracken_report", required=True)
    parser.add_argument("--target_rank", help="Target rank from e.g. ['D','P','C','O','F','G','S']", default="S")
    parser.add_argument("--min_num_reads", type=int, help="Threshold for min number reads", default=10)
    parser.add_argument("--max_n", type=int, help="Maximum number of tax ids to return", default=None)
    parser.add_argument("--min_percent", type=float, help="Threshold for min percentage of reads", default=1)

    args = parser.parse_args()

    tax_ids = get_top_n_kraken_hits(args.bracken_report,args.target_rank,args.min_num_reads,args.max_n, args.min_percent)
    for id in tax_ids:
        print(id)

if __name__ == "__main__":
    main()