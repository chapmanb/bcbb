"""Explore true/false WHAM calls for filtering approaches.
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import vcf

from bcbio import utils

def main(vcf_file, good_bed, bad_bed):
    class_csv = _convert_to_csv(vcf_file, good_bed, bad_bed)
    _plot_metrics(class_csv)

def _plot_metrics(class_csv):
    out_file = "%s.png" % os.path.splitext(class_csv)[0]
    sns.set(style="darkgrid")
    df = pd.read_csv(class_csv)
    g = sns.FacetGrid(df, row="class", col="attr", sharex=False, sharey=False, margin_titles=True)
    g.map(plt.hist, "val", color="steelblue", lw=0)
    g.fig.savefig(out_file)

def _convert_to_csv(vcf_file, good_bed, bad_bed):
    """Convert WHAM output file into BED format for graphical exploration.
    """
    attrs = ["PU", "LRT", "SI", "MQ"]
    buffer_size = 25  # bp around break ends
    out_file = "%s-metrics.csv" % utils.splitext_plus(vcf_file)[0]
    if not utils.file_uptodate(out_file, vcf_file):
        good = _read_bed(good_bed)
        bad = _read_bed(bad_bed)
        with open(out_file, "w") as out_handle:
            reader = vcf.Reader(filename=vcf_file)
            writer = csv.writer(out_handle)
            header = ["chrom", "start", "end", "class", "attr", "val"]
            writer.writerow(header)
            for rec in reader:
                start = max(rec.start - buffer_size, 0)
                if rec.INFO["BE"][0] not in [".", None]:
                    other_chrom, end, count = rec.INFO["BE"]
                    if int(end) > start and other_chrom == rec.CHROM:
                        end = int(end) + buffer_size
                        if (rec.CHROM, start, end) in good:
                            cur_class = "good"
                        elif (rec.CHROM, start, end) in bad:
                            cur_class = "bad"
                        else:
                            cur_class = None
                        if cur_class:
                            for attr in attrs:
                                writer.writerow([rec.CHROM, start, end, cur_class, attr, rec.INFO[attr]])
    return out_file

def _read_bed(in_file):
    out = set([])
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end = line.split("\t")[:3]
            out.add((chrom, int(start), int(end)))
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
