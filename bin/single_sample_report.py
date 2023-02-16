#!/usr/bin/env python
"""Create workflow report."""
# based of wf-metagenomics
import argparse
import json
import os

from aplanat import bars
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import layout
from bokeh.plotting import figure
from bokeh.resources import INLINE
from jinja2 import Template


def plot_hist_data(counts, edges, **kwargs):
    """Create histogram from json."""
    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(**kwargs)
    p = figure(**defaults)
    p.quad(
            top=counts, bottom=0, left=edges[:-1], right=edges[1:],
            alpha=0.6)
    return p


def kraken(summaries, section):
    """Kraken quality report section."""
    with open(summaries) as f:
        datas = json.load(f)
    bc_counts = {}
    for sample_id, data in datas.items():
        len_hist = list(data["len"].items())
        lbins, lcounts = list(zip(*len_hist))
        len_plot = plot_hist_data(
            lcounts, lbins, title="Read length distribution.",
            x_axis_label='Read Length / bases',
            y_axis_label='Number of reads'
        )
        qual_hist = list(data["qual"].items())
        edges, counts = list(zip(*qual_hist))
        total_reads = data["total_reads"]
        qual_plot = plot_hist_data(
            counts, edges, title="Read quality score",
            x_axis_label="Quality score",
            y_axis_label="Number of reads")

        section.markdown("##Read stats")
        section.plot(
            layout([[len_plot, qual_plot]], sizing_mode="stretch_width"))
        bc_counts[sample_id] = total_reads


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", required=True,
        help="Read summary JSON file.")
    parser.add_argument(
        "--lineages", nargs='+', required=True,
        help="Read lineage file.")
    parser.add_argument(
        "--vistempl", required=True)
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    args = parser.parse_args()

    report = WFReport(
        "Scylla Metagenomics Report for ", "scylla",
        revision=args.revision, commit=args.commit)

    templ = None
    with open(args.vistempl, "r") as vistempl:
        templ = vistempl.read()

    #
    # Sankey plot
    #
    all_json = {}
    for json_file in args.lineages:
        with open(json_file, 'r') as json_handle:
            sample_lineages = json.load(json_handle)
            all_json.update(sample_lineages)

    templ = templ.replace(
        "replace_me_data",
        json.dumps(all_json).replace('"', '\\"'))


    report.template = Template(templ)
    bokeh_resources = INLINE.render()
    report.template.render(bokeh_resources=bokeh_resources)

    #
    # Standard read metrics
    #
    section = report.add_section()
    kraken(args.stats, section)

    report.write(args.report)


if __name__ == "__main__":
    main()
