#!/usr/bin/env python
"""Convert a GFF2 file to an updated GFF3 format file.

Usage:
    gff2_to_gff3.py <in_gff2_file>

The output file has the same name with the extension gff3.
"""
import sys
import os

from BCBio.GFF.GFFParser import GFFAddingIterator
from BCBio.GFF.GFFOutput import GFF3Writer

def main(in_file):
    base, ext = os.path.splitext(in_file)
    out_file = "%s.gff3" % (base)
    in_handle = open(in_file)
    out_handle = open(out_file, "w")
    reader = GFFAddingIterator()
    writer = GFF3Writer()
    for recs in reader.get_features(in_handle, target_lines = 10000):
        writer.write(recs.values(), out_handle)

    in_handle.close()
    out_handle.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print __doc__
        sys.exit()
    main(sys.argv[1])
