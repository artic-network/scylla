#!/usr/bin/env python

import csv
from pathlib import Path
from mako.lookup import TemplateLookup
from datetime import date
import argparse
from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import json
import os


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

def get_binned_data(list_vals, num_bins, start=0):
    max_val = max(list_vals) + 1
    print (max_val)
    step = int((max_val - start) / num_bins) + ((max_val - start) % num_bins > 0)
    if step > 100:
        step = (int(step/100) + 1) *100
    print(step)
    binned_data = []
    for i in range(num_bins):
        binned_data.append({'bin_start':start + step*i, 'bin_end': start + step*(i+1), 'count':0})
    for j in list_vals:
        for bin in binned_data:
            if bin['bin_start'] <= j < bin['bin_end']:
                bin['count'] += 1
    print(binned_data)
    return binned_data, step

def summarize_read_counts(read_counts_file, num_bins=42, start=0):
    with open(read_counts_file.resolve(), "rt") as qc_file:
        reader = csv.DictReader(qc_file, delimiter="\t")
        lens = []
        quals = []
        for row in reader:
            lens.append(int(row["read_length"]))
            quals.append(float(row["mean_quality"]))
    len_data, len_step =  get_binned_data(lens, num_bins, start)
    qual_data, qual_step = get_binned_data(quals, num_bins, start)
    return  len_data, len_step, qual_data 


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--assignments",
        help="JSON file(s) of kraken/bracken assignments",
        nargs="+",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--read_counts", help="JSON file of read_counts", required=False, type=Path
    )
    parser.add_argument(
        "--warnings", help="text file containing any warnings", required=False, nargs='*', type=Path
    )
    parser.add_argument(
        "--sample_id", help="Unique ID of sample", required=False, default="sample"
    )

    parser.add_argument("--prefix", help="HTML output prefix ", default="scylla")

    parser.add_argument(
        "--template", help="HTML template for report", default="scylla.mako.html"
    )
    parser.add_argument("--classifier", help="Classifier used", default="Kraken")
    parser.add_argument("--classification_database", help="Database used for classification", default="PlusPF")
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
                assignments = assignments[:-1] + ", " + contents[1:]

    if args.read_counts:
        read_length_counts, read_length_step, read_quality_counts = summarize_read_counts(args.read_counts)
    else:
        read_length_counts, read_length_step, read_quality_counts = [], 2, []

    warnings = ""
    if len(args.warnings) > 0:
        for warning_file in args.warnings:
            if os.path.getsize(warning_file) > 0:
                with open(warning_file, "r") as f:
                    content = json.load(f)
                    warnings += content["msg"]

    data_for_report = {"sankey_data": assignments, "read_length_data": read_length_counts, "read_length_step": read_length_step, "read_quality_data": read_quality_counts, "classifier": args.classifier, "classification_database": args.classification_database, "warnings": warnings}

    outfile = args.prefix + "_report.html"
    make_output_report(outfile, args.template, args.version, sample, data_for_report)


if __name__ == "__main__":
    main()
