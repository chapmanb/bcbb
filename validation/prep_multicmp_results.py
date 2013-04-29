#!/usr/bin/env python
"""Prepare subsets of data from multi-method variant calling comparison.

Usage:
  prep_multicmp_results.py <in_csv>
"""
import bisect
import os
import sys

import pandas as pd

def main(in_file):
    out_file = "%s-prep%s" % os.path.splitext(in_file)
    df = pd.read_csv(in_file)
    df["aligner"] = [get_aligner(x) for x in df["sample"]]
    df["bamprep"] = [get_bamprep(x) for x in df["sample"]]
    df["sample"] = ["%s %s" % (x, y) if x and y else s
                    for (x, y, s) in zip(df["aligner"], df["bamprep"], df["sample"])]
    floors = get_group_floors(df)
    df["value.floor"] = [get_floor_value(x, cat, vartype, floors)
                         for (x, cat, vartype) in zip(df["value"], df["category"], df["variant.type"])]
    df.to_csv(out_file)

def get_floor_value(x, cat, vartype, floors):
    base = floors[(cat, vartype)]
    print cat, vartype, x, base
    return x - base

def get_group_floors(df):
    """Floor values to nearest 10,000 for each category.
    """
    floors = {}
    floor_vals = [x * 1e4 for x in range(50)]
    for name, group in df.groupby(["category", "variant.type"]):
        floors[name] = int(floor_vals[bisect.bisect(floor_vals, min(group["value"])) - 1])
    return floors

def get_aligner(x):
    if x in ["NA12878-1", "NA12878-3"]:
        return "novoalign"
    elif x in ["NA12878-2", "NA12878-4"]:
        return "bwa"
    else:
        return ""

def get_bamprep(x):
    if x in ["NA12878-1", "NA12878-2"]:
        return "gatk"
    elif x in ["NA12878-3", "NA12878-4"]:
        return "gkno"
    else:
        return ""

if __name__ == "__main__":
    main(sys.argv[1])
