#!/usr/bin/env python
"""Plot depth for a set of heterozygous calls relative to quality and allele ratio.

Used to help identify cutoff for filtering false positives by comparing
distribution to true positives.

Usage:
  plot_depth_ratio.py <VCF file of het calls> '<Plot title>'
"""
import os
import sys

import vcf
import prettyplotlib as ppl
import matplotlib.pyplot as plt

def main(in_file, title):
    depths, ratios, quals = get_ad_depth(in_file)
    plot_qual_hist(quals, in_file)
    plot_depth_ratios(depths, ratios, quals, in_file, title)

def plot_depth_ratios(depths, ratios, quals, in_file, title):
    out_file = "%s-depthratios.png" % os.path.splitext(in_file)[0]
    fig, ax = plt.subplots(1)
    for ds, rs, qualrange in _group_ratios_by_qual(depths, ratios, quals):
        print qualrange, len(ds)
        ppl.scatter(ax, x=depths, y=ratios, label=qualrange)
    ppl.legend(ax, title="Quality score range")
    ax.set_title(title)
    ax.set_xlabel("Depth")
    ax.set_ylabel("Variant/Total ratio")
    fig.savefig(out_file)

def _group_ratios_by_qual(depths, ratios, quals):
    #ranges = [(0, 100), (100, 250), (250, 500), (500, 1000), (1000, 2500)]
    #ranges = [(0, 50), (50, 100), (100, 150), (150, 250)]
    ranges = [(0, 250), (250, 500)]
    for qs, qe in ranges:
        cur_ds = []
        cur_rs = []
        for d, r, q in zip(depths, ratios, quals):
            if q >= qs and q < qe:
                cur_ds.append(d)
                cur_rs.append(r)
        yield cur_ds, cur_rs, "%s-%s" % (qs, qe)

def plot_qual_hist(quals, in_file):
    quals = [x for x in quals if x < 500.0]
    out_file = "%s-hist.png" % os.path.splitext(in_file)[0]
    fig, ax = plt.subplots(1)
    ppl.hist(ax, [quals], bins=100)
    fig.savefig(out_file)

def get_ad_depth(in_file):
    depths = []
    ratios = []
    quals = []
    with open(in_file) as in_handle:
        reader = vcf.Reader(in_handle)
        for rec in reader:
            for sample in rec.samples:
                try:
                    ad = sample["AD"]
                except AttributeError:
                    ad = []
                if len(ad) == 2:
                    ref, alt = sample["AD"]
                    depth = ref + alt
                    if depth > 0:
                        depths.append(min(rec.INFO["DP"], 500))
                        ratios.append(alt / float(depth))
                        quals.append(rec.QUAL)
    return depths, ratios, quals

if __name__ == "__main__":
    main(*sys.argv[1:])
