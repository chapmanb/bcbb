"""Test decoration of existing SeqRecords with GFF through a SeqIO interface.
"""
import sys
import os
import unittest
import pprint

from Bio import SeqIO
from BCBio.GFF.GFFParser import GFFMapReduceFeatureAdder

class MapReduceGFFTest(unittest.TestCase):
    """Tests GFF parsing using a map-reduce framework for parallelization.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.getcwd(), "GFF")
        self._test_gff_file = os.path.join(self._test_dir,
                "c_elegans_WS199_shortened_gff.txt")
        self._disco_host = "http://localhost:7000"
    
    def t_local_map_reduce(self):
        """General map reduce framework without parallelization.
        """
        cds_limit_info = dict(
                gff_types = [('Non_coding_transcript', 'gene'),
                             ('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        feature_adder = GFFMapReduceFeatureAdder(dict(), None)
        feature_adder.add_features(self._test_gff_file, cds_limit_info)
        final_rec = feature_adder.base['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature

    def t_disco_map_reduce(self):
        """Map reduce framework parallelized using disco.
        """
        # this needs to be more generalized but fails okay with no disco
        try:
            import disco
            import simplejson
        except ImportError:
            print "Skipping -- disco and json not found"
            return
        cds_limit_info = dict(
                gff_types = [('Non_coding_transcript', 'gene'),
                             ('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        feature_adder = GFFMapReduceFeatureAdder(dict(), self._disco_host)
        feature_adder.add_features(self._test_gff_file, cds_limit_info)
        final_rec = feature_adder.base['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature

class CElegansGFFTest(unittest.TestCase):
    """Real life test case using C elegans chromosome and GFF data

    Uses GFF3 data from:

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
        self._test_gff_ann_file = os.path.join(self._test_dir,
                "c_elegans_WS199_ann_gff.txt")
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
        seq_handle = open(seq_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        seq_handle.close()
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
        feature_adder.add_features(gff_file, limit_info)

    def _get_feature_adder(self):
        """Internal reusable function to get the feature adder.
        """
        seq_handle = open(self._test_seq_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        seq_handle.close()
        return GFFMapReduceFeatureAdder(seq_dict)
    
    def t_possible_limits(self):
        """Calculate possible queries to limit a GFF file.
        """
        feature_adder = self._get_feature_adder()
        possible_limits = feature_adder.available_limits(self._test_gff_file)
        print
        pprint.pprint(possible_limits)

    def t_flat_features(self):
        """Check addition of flat non-nested features to multiple records.
        """
        feature_adder = self._get_feature_adder()
        pcr_limit_info = dict(
            gff_types = [('Orfeome', 'PCR_product'),
                         ('GenePair_STS', 'PCR_product'),
                         ('Promoterome', 'PCR_product')]
            )
        feature_adder.add_features(self._test_gff_file, pcr_limit_info)
        assert len(feature_adder.base['I'].features) == 4
        assert len(feature_adder.base['X'].features) == 5

    def t_nested_features(self):
        """Check three-deep nesting of features with gene, mRNA and CDS.
        """
        feature_adder = self._get_feature_adder()
        cds_limit_info = dict(
                gff_types = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        feature_adder.add_features(self._test_gff_file, cds_limit_info)
        final_rec = feature_adder.base['I']
        # first gene feature is plain
        assert len(final_rec.features) == 2 # two gene feature
        assert len(final_rec.features[0].sub_features) == 1 # one transcript
        # 15 final CDS regions
        assert len(final_rec.features[0].sub_features[0].sub_features) == 15

    def t_nested_multiparent_features(self):
        """Verify correct nesting of features with multiple parents.
        """
        feature_adder = self._get_feature_adder()
        cds_limit_info = dict(
                gff_types = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        feature_adder.add_features(self._test_gff_file, cds_limit_info)
        final_rec = feature_adder.base['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature
        cur_subs = final_rec.features[1].sub_features
        assert len(cur_subs) == 3 # three transcripts
        # the first and second transcript have the same CDSs
        assert len(cur_subs[0].sub_features) == 6
        assert len(cur_subs[1].sub_features) == 6
        assert cur_subs[0].sub_features[0] is cur_subs[1].sub_features[0]

    def t_no_dict_error(self):
        """Ensure an error is raised when no dictionary to map to is present.
        """
        feature_adder = GFFMapReduceFeatureAdder(dict(), create_missing=False)
        try:
            feature_adder.add_features(self._test_gff_file)
            # no error -- problem
            raise AssertionError('Did not complain with missing dictionary')
        except KeyError:
            pass

    def t_gff_annotations(self):
        """Check GFF annotations placed on an entire sequence.
        """
        feature_adder = GFFMapReduceFeatureAdder(dict())
        feature_adder.add_features(self._test_gff_ann_file)
        final_rec = feature_adder.base['I']
        print final_rec.annotations

class SolidGFFTester(unittest.TestCase):
    """Test reading output from SOLiD analysis, as GFF3.

    See more details on SOLiD GFF here:

    http://solidsoftwaretools.com/gf/project/matogff/
    """
    def setUp(self):
        self._test_dir = os.path.join(os.getcwd(), "GFF")
        self._test_gff_file = os.path.join(self._test_dir,
                "F3-unique-3.v2.gff")

    def t_basic_solid_parse(self):
        """Basic parsing of SOLiD GFF results files.
        """
        feature_adder = GFFMapReduceFeatureAdder(dict())
        feature_adder.add_features(self._test_gff_file)
        test_feature = feature_adder.base['3_341_424_F3'].features[0]
        assert test_feature.location.nofuzzy_start == 102716
        assert test_feature.location.nofuzzy_end == 102736
        assert test_feature.qualifiers['score'] == ['10.6']
        assert test_feature.qualifiers['source'] == ['solid']
        assert test_feature.strand == -1
        assert test_feature.type == 'read'
        assert test_feature.qualifiers['g'] == ['T2203031313223113212']
        assert len(test_feature.qualifiers['q']) == 20

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
    tests = [CElegansGFFTest, MapReduceGFFTest, SolidGFFTester]
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)
    return test_suite

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
