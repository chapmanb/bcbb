"""Test decoration of existing SeqRecords with GFF through a SeqIO interface.
"""
from __future__ import with_statement
import sys
import os
import unittest
import pprint

from Bio import SeqIO
from BCBio.SeqIO.GFFIO import GFFFeatureAdder

class CElegansGFFTest(unittest.TestCase):
    """Real life test case using C elegans chromosome and GFF data

    Uses data from:

    ftp://ftp.wormbase.org/pub/wormbase/genomes/c_elegans/
    genome_feature_tables/GFF3/
    ftp://ftp.wormbase.org/pub/wormbase/genomes/c_elegans/sequences/dna/
    """
    def setUp(self):
        self._test_dir = os.path.join(os.getcwd(), "GFF")
        self._test_seq_file = os.path.join(self._test_dir,
                "c_elegans_WS199_dna_shortened.fa")
        self._test_gff_file = os.path.join(self._test_dir,
                "c_elegans_WS199_shortened_gff.txt")
        self._full_dir = "/usr/home/chapmanb/mgh/ruvkun_rnai/wormbase/" + \
                "data_files_WS198"

    def not_t_full_celegans(self):
        """Test the full C elegans chromosome and GFF files.

        This is used to test GFF on large files and is not run as a standard
        test. You will need to download the files and adjust the paths
        to run this.
        """
        # read the sequence information
        seq_file = os.path.join(self._full_dir, "c_elegans.WS199.dna.fa")
        gff_file = os.path.join(self._full_dir, "c_elegans.WS199.gff3")
        with open(seq_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        feature_adder = GFFFeatureAdder(seq_dict)
        #with open(gff_file) as gff_handle:
        #    possible_limits = feature_adder.available_limits(gff_handle)
        #    pprint.pprint(possible_limits)
        rnai_types = [('Orfeome', 'PCR_product'),
                    ('GenePair_STS', 'PCR_product'),
                    ('Promoterome', 'PCR_product')]
        gene_types = [('Non_coding_transcript', 'gene'),
                      ('Coding_transcript', 'gene'),
                      ('Coding_transcript', 'mRNA'),
                      ('Coding_transcript', 'CDS')]
        limit_info = dict(gff_types = rnai_types + gene_types)
        with open(gff_file) as gff_handle:
            feature_adder.add_features(gff_handle, limit_info)

    def t_possible_limits(self):
        """Calculate possible queries to limit a GFF file.
        """
        with open(self._test_seq_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        feature_adder = GFFFeatureAdder(seq_dict)
        with open(self._test_gff_file) as gff_handle:
            possible_limits = feature_adder.available_limits(gff_handle)
            print
            pprint.pprint(possible_limits)

    def t_flat_features(self):
        """Check addition of flat non-nested features to multiple records.
        """
        with open(self._test_seq_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        feature_adder = GFFFeatureAdder(seq_dict)
        pcr_limit_info = dict(
            gff_types = [('Orfeome', 'PCR_product'),
                         ('GenePair_STS', 'PCR_product'),
                         ('Promoterome', 'PCR_product')]
            )
        with open(self._test_gff_file) as gff_handle:
            feature_adder.add_features(gff_handle, pcr_limit_info)
        assert len(feature_adder.base['I'].features) == 4
        assert len(feature_adder.base['X'].features) == 5

    def t_nested_features(self):
        """Check three-deep nesting of features with gene, mRNA and CDS.
        """
        with open(self._test_seq_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        feature_adder = GFFFeatureAdder(seq_dict)
        cds_limit_info = dict(
                gff_types = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        with open(self._test_gff_file) as gff_handle:
            feature_adder.add_features(gff_handle, cds_limit_info)
        final_rec = feature_adder.base['I']
        # first gene feature is plain
        assert len(final_rec.features) == 2 # two gene feature
        assert len(final_rec.features[0].sub_features) == 1 # one transcript
        # 15 final CDS regions
        assert len(final_rec.features[0].sub_features[0].sub_features) == 15

    def t_nested_multiparent_features(self):
        """Verify correct nesting of features with multiple parents.
        """
        with open(self._test_seq_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        feature_adder = GFFFeatureAdder(seq_dict)
        cds_limit_info = dict(
                gff_types = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        with open(self._test_gff_file) as gff_handle:
            feature_adder.add_features(gff_handle, cds_limit_info)
        final_rec = feature_adder.base['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature
        cur_subs = final_rec.features[1].sub_features
        assert len(cur_subs) == 3 # three transcripts
        # the first and second transcript have the same CDSs
        assert len(cur_subs[0].sub_features) == 6
        assert len(cur_subs[1].sub_features) == 6
        assert cur_subs[0].sub_features[0] is cur_subs[1].sub_features[0]

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
    tests = [CElegansGFFTest]
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)
    return test_suite

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
