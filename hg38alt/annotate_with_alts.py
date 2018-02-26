#!/usr/bin/env python
"""Annotate BED file with alts overlapping regions.

Usage:
  annotate_with_alts.py <orig_bed> <alt_bed>
"""
from __future__ import print_function
import os
import sys

import pybedtools

def main(orig_bed, alt_bed):
    last = None
    cur_alt = []
    out_file = "%s-withalts.bed" % (os.path.splitext(orig_bed)[0])
    stats = CollectStats()
    with open(out_file, "w") as out_handle:
        for r in pybedtools.BedTool(orig_bed).intersect(alt_bed, wao=True):
            orig = tuple(r.fields[:-5])
            alt = r.fields[-2] if r.fields[-2] != "." else None
            if alt and orig == last:
                cur_alt.append(alt)
            else:
                if last:
                    alts = ",".join([a for a in cur_alt if a])
                    stats.add(orig, alts)
                    if not alts:
                        alts = "."
                    out_handle.write("\t".join(list(last) + [alts]) + "\n")
                last = orig
                cur_alt = [alt]
        if last:
            alts = ",".join([a for a in cur_alt if a])
            stats.add(last, alts)
            if not alts:
                alts = "."
            out_handle.write("\t".join(list(last) + [alts]) + "\n")
    stats.print()

class CollectStats:
    def __init__(self):
        self.genes = set([])
        self.total_regions = 0
        self.total_bases = 0
        self.alt_bases = 0
        self.alt_regions = 0

    def add(self, region, alt):
        cur_size = int(region[2]) - int(region[1])
        self.total_bases += cur_size
        self.total_regions += 1
        if alt:
            self.alt_regions += 1
            self.alt_bases += cur_size
            self.genes.add(region[3])

    def print(self):
        print("With alts: %s out of %s regions %0.1f%%; %sbp out of %sbp: %0.1f%%" %
              (self.alt_regions, self.total_regions,
               float(self.alt_regions) / self.total_regions * 100.0,
               self.alt_bases, self.total_bases,
               float(self.alt_bases) / self.total_bases * 100.0))
        for g in  sorted(list(self.genes)):
            print(g)

if __name__ == "__main__":
    main(*sys.argv[1:])
