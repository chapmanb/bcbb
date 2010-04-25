#!/usr/bin/env python
"""Sort MAF file by the size of alignments -- largest to smallest.

Usage:
    maf_sort_by_size.py <maf file>
"""
from __future__ import with_statement
import sys
import os

from bx.align import maf
from bx import interval_index_file

def main(in_file):
    base, ext = os.path.splitext(in_file)
    out_file = "%s-sorted%s" % (base, ext)
    index_file = in_file + ".index"
    if not os.path.exists(index_file):
        build_index(in_file, index_file)

    # pull out the sizes and positions of each record
    rec_info = []
    with open(in_file) as in_handle:
        reader = maf.Reader(in_handle)
        while 1:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            rec_info.append((rec.text_size, pos))
    rec_info.sort(reverse=True)

    # write the records in order, pulling from the index
    index = maf.Indexed(in_file, index_file)
    with open(out_file, "w") as out_handle:
        writer = maf.Writer(out_handle)
        for size, pos in rec_info:
            rec = index.get_at_offset(pos)
            writer.write(rec)

def build_index(in_file, index_file):
    """Build an index of the MAF file for retrieval.
    """
    indexes = interval_index_file.Indexes()
    with open(in_file) as in_handle:
        reader = maf.Reader(in_handle)
        while 1:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            for c in rec.components:
                indexes.add(c.src, c.forward_strand_start,
                        c.forward_strand_end, pos, max=c.src_size )

    with open(index_file, "w") as index_handle:
        indexes.write(index_handle)

if __name__ == "__main__":
    main(sys.argv[1])
