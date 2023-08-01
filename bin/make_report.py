#!/usr/bin/env python

import json
import csv
from collections import defaultdict
from pathlib import Path
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date
import argparse
from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import pandas as pd
import os
from Bio import SeqIO


def make_output_report(
    report_to_generate, template, version, sample, data_for_report={}
):
    template_dir = os.path.abspath(os.path.dirname(__file__))
    mylookup = TemplateLookup(directories=[template_dir])  # absolute or relative works
    mytemplate = mylookup.get_template(template)

    buf = StringIO()

    print(template_dir)
    js_path = os.path.join(template_dir, "report_utils", "sankey.js")
    print(js_path)
    with open(js_path, "r") as f:
        data_for_report["sankey_js"] = f.read()

    ctx = Context(
        buf,
        date=date.today(),
        version=version,
        sample=sample,
        data_for_report=data_for_report,
    )

    try:
        mytemplate.render_context(ctx)

    except:
        traceback = RichTraceback()
        for filename, lineno, function, line in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, "w") as fw:
        print("Generating: " + f"{report_to_generate}")
        fw.write(buf.getvalue())


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--assignments",
        help="JSON file(s) of kraken/bracken assignments",
        nargs='+',
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--read_counts", help="JSON file of read_counts", required=False, type=Path
    )
    parser.add_argument(
        "--sample_id", help="Unique ID of sample", required=False, default="sample"
    )

    parser.add_argument("--prefix", help="HTML output prefix ", default="scylla")

    parser.add_argument(
        "--template", help="HTML template for report", default="scylla.mako.html"
    )
    parser.add_argument("--version", help="Scylla version", default="unknown")

    args = parser.parse_args()

    sample = args.sample_id

    assignments = None
    for bracken_file in args.assignments:
        with open(bracken_file.resolve(), "rt") as bracken_handle:
            contents = bracken_handle.read().strip().replace('"', '\\"')
            if not assignments:
                assignments = contents
            else:
                assignments = assignments[:-1] + ", " + contents

    read_counts = ""
    if args.read_counts:
        with open(read_counts.resolve(), "rt") as counts_handle:
            read_counts = counts_handle.read().strip().replace('"', '\\"')
    data_for_report = {"sankey_data": assignments, "read_count_data": read_counts}

    outfile = args.prefix + "_report.html"
    make_output_report(outfile, args.template, args.version, sample, data_for_report)


if __name__ == "__main__":
    main()
