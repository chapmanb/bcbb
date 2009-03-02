#!/usr/bin/env python
"""Examine conservation of a protein by comparison to BLAST hits.

Given a UniProt protein ID or accession number as input (really anything that
can be queried in NCBI), this performs a BLAST search against the
non-redundant protein database and parses the results. Using them, a plot is
generated of average conservation across the protein. This provides a quick
evaluation of conserved and fluctuating regions.

Usage:
    blast_conservation_plot.py <accession>
"""
from __future__ import with_statement
import sys
import os

from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SubsMat import MatrixInfo
import pylab
import numpy

def main(accession):
    window_size = 29
    cache_dir = os.path.join(os.getcwd(), "cache")
    ncbi_manager = NCBIManager(cache_dir)
    protein_gi = ncbi_manager.search_for_gi(accession, "protein")
    blast_rec = ncbi_manager.remote_blast(protein_gi, "blastp")
    cons_caculator = BlastConservationCalculator()
    data_smoother = SavitzkyGolayDataSmoother(window_size)
    cons_dict = cons_caculator.conservation_dict(blast_rec)
    indexes = cons_dict.keys()
    indexes.sort()
    pos_data = []
    cons_data = []
    for pos in indexes:
        pos_data.append(pos + 1)
        if len(cons_dict[pos]) > 0:
            cons_data.append(numpy.median(cons_dict[pos]))
        else:
            cons_data.append(0)
    smooth_data = data_smoother.smooth_values(cons_data)
    smooth_pos_data = pos_data[data_smoother.half_window():
            len(pos_data) - data_smoother.half_window()]
    pylab.plot(smooth_pos_data, smooth_data)
    pylab.axis(xmin=min(pos_data), xmax=max(pos_data))
    pylab.xlabel("Amino acid position")
    pylab.ylabel("Conservation")
    pylab.savefig('%s_conservation.png' % accession.replace(".", "_"))

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

class BlastConservationCalculator:
    """Calculate conservation across a protein from a BLAST record.
    """
    def __init__(self, matrix_name="blosum62"):
        """Initialize with the name of a substitution matrix for comparisons.
        """
        self._subs_mat = getattr(MatrixInfo, matrix_name)
        self._no_use_thresh = 0.95

    def conservation_dict(self, blast_rec):
        """Get dictionary containing substitution scores based on BLAST HSPs.
        """
        cons_dict = {}
        rec_size = int(blast_rec.query_letters)
        for base_index in range(rec_size):
            cons_dict[base_index] = []
        for align in blast_rec.alignments:
            for hsp in align.hsps:
                if (float(hsp.identities) / float(rec_size) <=
                        self._no_use_thresh):
                    cons_dict = self._add_hsp_conservation(hsp, cons_dict)
        return cons_dict

    def _add_hsp_conservation(self, hsp, cons_dict):
        """Add conservation information from an HSP BLAST alignment.
        """
        start_index = int(hsp.query_start) - 1
        hsp_index = 0
        for q_index in range(len(hsp.query)):
            if (hsp.query[q_index] != '-'):
                if (hsp.sbjct[q_index] != '-'):
                    try:
                        sub_val = self._subs_mat[(hsp.query[q_index],
                                                  hsp.sbjct[q_index])]
                    except KeyError:
                        sub_val = self._subs_mat[(hsp.sbjct[q_index],
                                                  hsp.query[q_index])]
                    cons_dict[start_index + hsp_index].append(sub_val)
                hsp_index += 1
        return cons_dict

class NCBIManager:
    """Manage interactions with NCBI through Biopython
    """
    def __init__(self, cache_dir):
        self._cache_dir = cache_dir
        if not(os.path.exists(cache_dir)):
            os.makedirs(cache_dir)

    def search_for_gi(self, uniprot_id, db_name):
        """Find the NCBI GI number corresponding to the given input ID.
        """
        handle = Entrez.esearch(db=db_name, term=uniprot_id)
        record = Entrez.read(handle)
        ids = record["IdList"]
        if len(ids) == 0:
            raise ValueError("Not found in NCBI: %s" % ids)
        return ids[0]

    def remote_blast(self, search_gi, blast_method):
        """Perform a BLAST against the NCBI server, returning the record.
        """
        out_file = os.path.join(self._cache_dir, "%s_%s_blo.xml" % (blast_method,
            search_gi))
        if not os.path.exists(out_file):
            blast_handle = NCBIWWW.qblast(blast_method, "nr", search_gi)
            with open(out_file, 'w') as out_handle:
                for line in blast_handle:
                    out_handle.write(line)
            blast_handle.close()
        with open(out_file) as in_handle:
            rec_it = NCBIXML.parse(in_handle)
            return rec_it.next()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1])
