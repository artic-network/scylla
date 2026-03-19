#!/usr/bin/env python
import argparse
import gzip
import json
from pathlib import Path


def load_spike_in_dict(ref_file):
    with open(ref_file) as handle:
        return json.load(handle)


def concatenate_references(ref_paths, output_file):
    with open(output_file, "w") as f:
        for ref_path in ref_paths:
            with gzip.open(ref_path, "rt") as gz_file:
                for line in gz_file:
                    f.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Combine spike-in reference sequences"
    )
    parser.add_argument(
        "--spike_ins",
        required=True,
        help="Comma-separated list of spike-in names",
    )
    parser.add_argument(
        "--spike_in_dict",
        required=True,
        help="JSON file mapping spike-in names to reference files",
    )
    parser.add_argument(
        "--spike_in_ref_dir",
        required=True,
        help="Directory containing spike-in reference files",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="combined_spikes.fa",
        help="Output FASTA file",
    )

    args = parser.parse_args()

    spike_names = [name.strip() for name in args.spike_ins.split(",") if name.strip()]
    spike_map = load_spike_in_dict(args.spike_in_dict)
    base_dir = Path(args.spike_in_ref_dir)
    ref_paths = [str(base_dir / spike_map[name]["ref"]) for name in spike_names]
    concatenate_references(ref_paths, args.output)


if __name__ == "__main__":
    main()
