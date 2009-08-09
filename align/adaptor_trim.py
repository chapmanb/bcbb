#!/usr/bin/env python
"""Trim adaptor sequences from reads; designed for short read sequencing.

Allows trimming of adaptor sequences from a list of SeqRecords produced
by the Biopython SeqIO library.
"""
from decorator import decorator

from Bio import pairwise2
from Bio.Seq import Seq

def _remove_adaptor(seq, region, right_side=True):
    """Remove an adaptor region and all sequence to the right or left.
    """
    # A check for repetitive regions. We handle them below by searching
    # from the left or right depending on the trimming method which should
    # be the expected result.
    #size = len(region)
    #pieces = [str(seq[i:i+size]) for i in range(len(seq) - size)]
    #if pieces.count(region) != 1:
    #    raise ValueError("Non-single match: %s to %s" % (region, seq))
    if right_side:
        pos = seq.find(region)
        return seq[:pos]
    else:
        pos = seq.rfind(region)
        return seq[pos+len(region):]

@decorator
def accept_seq_or_str(fn, *args, **kws):
    """Decorator allowing either a string or Seq as the first argument.

    A string or Seq is returning, determined by the input argument. The
    decorated function works on a string.
    """
    # check our sequence passed in as a keyword or first argument
    seq = kws.get('seq', None) or args[0]
    if isinstance(seq, Seq):
        if kws.has_key('seq'):
            kws['seq'] = str(seq)
        else:
            args = (str(seq),) + args[1:]
    else:
        seq = None
    final_seq = fn(*args, **kws)
    if seq is not None:
        final_seq = Seq(final_seq, seq.alphabet)
    return final_seq

@accept_seq_or_str
def trim_adaptor(seq, adaptor, num_errors, right_side=True):
    gap_char = '-'
    seq_a, adaptor_a, score, start, end = pairwise2.align.localms(str(seq),
            str(adaptor), 5.0, -4.0, -10.0, -0.5,
            one_alignment_only=True, gap_char=gap_char)[0]
    #print seq_a, adaptor_a, score, start, end
    seq_region = seq_a[start:end]
    adapt_region = adaptor_a[start:end]
    diffs = sum((0 if s == adapt_region[i] else 1) for i, s in
            enumerate(seq_region))
    # too many errors -- no trimming
    if diffs > num_errors:
        return seq
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq, seq_region.replace(gap_char, ""),
                right_side)

# ------- Testing Code
import sys
import unittest

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

class AdaptorAlignTrimTest(unittest.TestCase):
    """Test remove adaptor sequences using local alignments.
    """
    def t_1_simple_trim(self):
        """Trim adaptor from non-complex region with errors and deletions.
        """
        adaptor = "GATCGATCGATC"
        tseq = trim_adaptor("GGG" + adaptor + "CCC", adaptor, 2) # exact
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GATCGTTCGATC" + "CCC", adaptor, 2)# 1 error
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GATCGTTCGAAC" + "CCC", adaptor, 2)# 2 errors
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GATCGATCGTC" + "CCC", adaptor, 2) # deletion
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GAACGTTGGATC" + "CCC", adaptor, 2)# 3 errors
        assert tseq == "GGGGAACGTTGGATCCCC"
        tseq = trim_adaptor("GGG" + "CATCGGACGTAT" + "CCC", adaptor, 2)# very bad
        assert tseq == "GGGCATCGGACGTATCCC"

    def t_2_alternative_side_trim(self):
        """Trim adaptor from both sides of the sequence.
        """
        adaptor = "GATCGATCGATC"
        tseq = trim_adaptor("GGG" + "GATCGTTCGATC" + "CCC", adaptor, 2)
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GATCGTTCGATC" + "CCC", adaptor, 2, False)
        assert tseq == "CCC"

    def t_3_repetitive_sequence(self):
        """Trim a repetitive adaptor sequence.
        """
        adaptor = "GATCGATC"
        tseq = trim_adaptor("GGG" + "GATCGATCGATC" + "CCC", adaptor, 2)
        assert tseq == "GGG"
        tseq = trim_adaptor("GGG" + "GATCGATCGATC" + "CCC", adaptor, 2, False)
        assert tseq == "CCC"

    def t_4_passing_seqs(self):
        """Handle both Biopython Seq objects and strings.
        """
        adaptor = "GATCGATCGATC"
        tseq = trim_adaptor(Seq("GGG" + "GATCGTTCGATC" + "CCC", unambiguous_dna),
            adaptor, 2)
        assert isinstance(tseq, Seq)
        assert tseq.alphabet == unambiguous_dna
        tseq = trim_adaptor(
                seq = Seq("GGG" + "GATCGTTCGATC" + "CCC", unambiguous_dna),
                adaptor=adaptor,
                num_errors=2)
        assert isinstance(tseq, Seq)
        assert tseq.alphabet == unambiguous_dna

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [AdaptorAlignTrimTest]
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)
    return test_suite

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
