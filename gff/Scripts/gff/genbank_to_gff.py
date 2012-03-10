#!/usr/bin/env python
"""Convert a GenBank file into GFF format.

Usage:
    genbank_to_gff.py <genbank_file>
"""
import sys
import os

from Bio import SeqIO
from Bio import Seq

from BCBio import GFF

def main(gb_file):
    out_file = "%s.gff" % os.path.splitext(gb_file)[0]
    with open(out_file, "w") as out_handle:
        GFF.write(SeqIO.parse(gb_file, "genbank"), out_handle)

if __name__ == "__main__":
    main(*sys.argv[1:])
