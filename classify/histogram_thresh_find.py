#!/usr/bin/env python
"""Classify proteins using visual inspection of a histogram.

This uses charge regions to classify, but could be modified to any
type of useful classification.

Usage:
    histogram_thresh_find.py <fasta positive records> <fasta negative records>
"""
from __future__ import with_statement
import sys

from Bio import SeqIO
from Bio.SeqUtils import IsoelectricPoint
from Bio.SeqUtils import ProtParam

import numpy
import pylab

def main(pos_file, neg_file):
    cur_window = 75
    pos_charges = file_charges(pos_file, cur_window)
    neg_charges = file_charges(neg_file, cur_window)
    print numpy.mean(pos_charges), numpy.mean(neg_charges)
    n, bins, patches = pylab.hist([pos_charges, neg_charges], 30,
            normed=True, histtype='bar')
    pylab.xlabel('Isoelectric point')
    pylab.ylabel('Normalized percent of regions')
    pylab.title('Protein charge of %s amino acid windows' % cur_window)
    pylab.legend([patches[0][0], patches[1][0]], ['positives', 'negatives'])
    pylab.savefig('pos_neg_hist.png')
    pylab.show()

def file_charges(in_file, cur_window):
    """Handle calculation of charges for all records in a file.
    """
    all_charges = []
    with open(in_file) as in_handle:
        for rec in SeqIO.parse(in_handle, "fasta"):
            cur_charges = calc_region_charges(rec.seq, cur_window)
            above_thresh = [c for c in cur_charges if c >= 10.2]
            if above_thresh:
                print rec.name, len(above_thresh) / float(len(rec.seq))
            all_charges.extend(cur_charges)
    return all_charges

def calc_region_charges(seq, cur_window):
    """Perform calculation of charges via isoelectric points for a sequence.
    """
    # internal small regions, so do not deal with C and N terminal charges
    IsoelectricPoint.pKcterminal = {}
    IsoelectricPoint.pKnterminal = {}
    cur_pos = 0
    region_charges = []
    while cur_pos < len(seq) - cur_window:
        cur_seq = seq[cur_pos:cur_pos + cur_window]
        prot_analysis = ProtParam.ProteinAnalysis(str(cur_seq))
        ie_calc = IsoelectricPoint.IsoelectricPoint(cur_seq,
                prot_analysis.count_amino_acids())
        region_charges.append(ie_calc.pi())
        cur_pos += 1
    return region_charges

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1], sys.argv[2])
