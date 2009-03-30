#!/usr/bin/env python
"""Combine predictions of intrinsically disordered regions in proteins.

There are a number of different methods to calculate intrinsically disordered
regions in proteins. A summary of methods from Wikipedia is here:

http://en.wikipedia.org/wiki/Intrinsically_unstructured_proteins

This parses results from a few of these and combines them into a final
probabilistic prediction of disorder. It handles output from the following
programs:

    Disembl -- http://dis.embl.de/
    GlobPlot -- http://globplot.embl.de/
    Disopred2 -- http://bioinf.cs.ucl.ac.uk/disopred/

Usage:
    intrinsic_disorder_combiner.py <input_prefix>

Where input_prefix specifies the base name of the input files. They should have
extensions according to the program, like:
    myproteins-disembl.txt
    myproteins-glob.txt
    myproteins-disopred.txt
You also need a protein sequence file named similarly:
    myproteins-seq.txt
"""
from __future__ import with_statement
import sys
import os
import glob
import collections
from optparse import OptionParser

from Bio import SeqIO
import pylab
import numpy

class ProteinIdPredictions:
    """Hold all predictions for individual bases from various programs.
    """
    def __init__(self):
        self.base_predictions = collections.defaultdict(list)

class GlobPlotParse:
    """Parser to handle the results from the GlobPlot standalone program.
    """
    def __init__(self):
        self._pname = 'GlobPlot'

    def add_predictions(self, cur_handle, id_preds):
        """Add GlobPlot predictions on disorder.
        """
        for rec in SeqIO.parse(cur_handle, "fasta"):
            uniprot_id, glob_info, disorder_info = rec.description.split('|')
            ranges = [r.split('-') for r in
                    disorder_info.split(':')[-1].split(",") if r]
            ranges = [(int(x) - 1, int(y)) for x, y in ranges]
            for s, e in ranges:
                for index in range(s, e + 1):
                    id_preds[uniprot_id].base_predictions[index].append(
                            self._pname)
        return id_preds

class DisemblParse:
    """Parser handling the results of Disembl analysis.
    """
    def __init__(self):
        pass
    
    def add_predictions(self, cur_handle, id_preds):
        for rec in SeqIO.parse(cur_handle, "fasta"):
            pieces = rec.description.split()
            uniprot_id, pname = pieces[0].split('_')
            ranges = [r.replace(',', '').split('-') for r in pieces[1:]]
            ranges = [(int(x) - 1, int(y)) for x, y in ranges]
            for s, e in ranges:
                for index in range(s, e + 1):
                    id_preds[uniprot_id].base_predictions[index].append(
                            pname)
        return id_preds

def plot_disorder(rec, base_predictions, smooth_window, ymax):
    preds = [len(base_predictions[i]) for i in range(len(rec.seq))]
    positions = [i + 1 for i in range(len(rec.seq))]
    data_smoother = SavitzkyGolayDataSmoother(smooth_window)
    smooth_data = data_smoother.smooth_values(preds)
    smooth_pos_data = positions[data_smoother.half_window():
            len(preds) - data_smoother.half_window()]
    pylab.clf()
    pylab.plot(smooth_pos_data, smooth_data)
    pylab.ylabel('Disorder support (number of predictions)')
    pylab.xlabel('Amino acid position')
    pylab.axis(ymin=0)
    pylab.ylim(ymin=0, ymax=ymax + .5)
    pylab.title('%s: Intrinsically disordered regions -- %s window' %
            (rec.id, smooth_window))
    file_name = "%s-idr.png"  % (rec.id)
    pylab.savefig(file_name)
    return file_name

def main(input_prefix):
    smooth_window = 75
    id_preds = collections.defaultdict(lambda : ProteinIdPredictions())
    program_info = [("glob", GlobPlotParse), ("disembl", DisemblParse)]
    for name_suffix, parse_class in program_info:
        parser = parse_class()
        cur_file = glob.glob("%s*%s*" % (input_prefix, name_suffix))[0]
        with open(cur_file) as cur_handle:
            id_preds = parser.add_predictions(cur_handle, id_preds)
    # find the highest number of predictions of all programs
    ymax = max(max([[len(preds) for preds in obj.base_predictions.values()] for
        obj in id_preds.values()]))
        
    seq_file = glob.glob("%s*seqs*" % (input_prefix))[0]
    with open(seq_file) as seq_handle:
        for rec in SeqIO.parse(seq_handle, "fasta"):
            plot_disorder(rec, id_preds[rec.id].base_predictions, smooth_window,
                    ymax)

class SavitzkyGolayDataSmoother:
    """Smooth data using the Savitzky-Golay technique from:

    http://www.dalkescientific.com/writings/NBN/plotting.html
    """
    def __init__(self, window_size):
        self._window_size = window_size
        if self._window_size%2 != 1:
            raise TypeError("smoothing requires an odd number of weights")

    def half_window(self):
        return (self._window_size-1)/2

    def smooth_values(self, values):
        half_window = (self._window_size-1)/2
        weights = self.savitzky_golay_weights(self._window_size)
        weights = [w*100.0 for w in weights]

        # Precompute the offset values for better performance.
        offsets = range(-half_window, half_window+1)
        offset_data = zip(offsets, weights)

        # normalize the weights in case the sum != 1
        total_weight = sum(weights)

        weighted_values = []
        for i in range(half_window, len(values)-half_window):
            weighted_value = 0.0
            for offset, weight in offset_data:
                weighted_value += weight*values[i+offset]
            weighted_values.append(weighted_value / total_weight)

        return weighted_values

    def savitzky_golay(self, window_size=None, order=2):
        if window_size is None:
            window_size = order + 2

        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window size is too small for the polynomial")

        # A second order polynomial has 3 coefficients
        order_range = range(order+1)
        half_window = (window_size-1)//2
        B = numpy.array(
            [ [k**i for i in order_range] for k in range(-half_window, half_window+1)] )

        #           -1
        # [  T     ]      T
        # [ B  * B ]  *  B
        M = numpy.dot(
               numpy.linalg.inv(numpy.dot(numpy.transpose(B), B)),
               numpy.transpose(B)
               )
        return M

    def savitzky_golay_weights(self, window_size=None, order=2, derivative=0):
        # The weights are in the first row
        # The weights for the 1st derivatives are in the second, etc.
        return self.savitzky_golay(window_size, order)[derivative]

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print __doc__
        sys.exit()
    main(args[0])
