"""Access an GFF file using bx-python's interval indexing.

Requires:
    bx-python: http://bitbucket.org/james_taylor/bx-python/wiki/Home
    gff library: http://github.com/chapmanb/bcbb/tree/master/gff

Index time:
  44 Mb file
  11 seconds
  Index is 7.5Mb
"""
from __future__ import with_statement
import os
import sys

from bx import interval_index_file

from BCBio import GFF

def main(gff_file):
    gff_index = gff_file + ".index"
    if not os.path.exists(gff_index):
        print "Indexing GFF file"
        index(gff_file)
    index = GFFIndexedAccess(gff_file, keep_open=True)
    print index.seqids
    print
    for feature in index.get_features_in_region("Chr2", 17500, 20000):
        print feature
    for feature in index.get_features_in_region("Chr5", 500000, 502500):
        print feature

    exam = GFF.GFFExaminer()
    #print exam.available_limits(gff_file)
    #print exam.parent_child_map(gff_file)

    found = 0
    limit_info = dict(
            gff_type = ["protein", "gene", "mRNA", "exon", "CDS", "five_prime_UTR",
                "three_prime_UTR"]
            )
    for feature in index.get_features_in_region("Chr1", 0, 50000, 
            limit_info):
        found += 1
    print found

class GFFIndexedAccess(interval_index_file.AbstractIndexedAccess):
    """Provide indexed access to a GFF file.
    """
    def __init__(self, *args, **kwargs):
        interval_index_file.AbstractIndexedAccess.__init__(self, *args,
                **kwargs)
        self._parser = GFF.GFFParser()

    @property
    def seqids(self):
        return self.indexes.indexes.keys()

    def get_features_in_region(self, seqid, start, end, limit_info=None):
        """Retrieve features located on a given region in start/end coordinates.
        """
        limit_info = self._parser._normalize_limit_info(limit_info)
        line_gen = self.get_as_iterator(seqid, int(start), int(end))
        recs = None
        for results in self._parser._lines_to_out_info(line_gen, limit_info):
            assert not recs, "Unexpected multiple results"
            recs = self._parser._results_to_features(dict(), results)
        if recs is None:
            return []
        else:
            assert len(recs) == 1
            rec = recs[seqid]
            return rec.features

    def read_at_current_offset(self, handle, **kwargs):
        line = handle.readline()
        return line

def index(gff_file, index_file=None):
    index = interval_index_file.Indexes()
    with open(gff_file) as in_handle:
        while 1:
            pos = in_handle.tell()
            line = in_handle.readline()
            if not line:
                break
            if not line.startswith("#"):
                parts = line.split("\t")
                (seqid, gtype, source, start, end) = parts[:5]
                index.add(seqid, int(start), int(end), pos)
    if index_file is None:
        index_file = gff_file + ".index"
    with open(index_file, "w") as index_handle:
        index.write(index_handle)
    return index_file

if __name__ == "__main__":
    main(*sys.argv[1:])
