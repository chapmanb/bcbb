"""Tests associated with detecting sequencing results dumped from a machine.
"""
import os
import unittest

import yaml

from bcbio.solexa import samplesheet

class SampleSheetTest(unittest.TestCase):
    """Deal with Illumina SampleSheets and convert to YAML input.
    """
    def setUp(self):
        self.ss_file = os.path.join(os.path.dirname(__file__),
                                    "data", "samplesheets",
                                    "illumina_samplesheet.csv")
        self.ss_non_multiplex = os.path.join(os.path.dirname(__file__),
                                             "data", "samplesheets",
                                             "illumina_samplesheet_non_multiplex_samples.csv")
        self.out_file = ""

    def tearDown(self):
        if os.path.exists(self.out_file):
            os.remove(self.out_file)

    def test_toyaml(self):
        """Convert CSV Illumina SampleSheet to YAML.
        """
        info = self.toyaml(self.ss_file)
        assert info[0]['lane'] == '1'
        assert info[0]['multiplex'][0]['barcode_id'] == 5

    def test_non_multiplexed(self):
        info = self.toyaml(self.ss_non_multiplex)
        assert os.path.exists(self.out_file)
        with open(self.out_file) as in_handle:
            info = yaml.load(in_handle)

        assert info[7]['lane'] == '8'
        assert not info[7].has_key("multiplex")

    def test_checkforrun(self):
        """Check for the presence of runs in an Illumina SampleSheet.
        """
        fcdir = "fake/101007_80HM7ABXX"
        config = {"samplesheet_directories" : [os.path.dirname(self.ss_file)]}
        ss = samplesheet.run_has_samplesheet(fcdir, config, False)
        assert ss is not None
        fcdir = "fake/101007_NOPEXX"
        ss = samplesheet.run_has_samplesheet(fcdir, config, False)
        assert ss is None


# Helper functions
# 

    def toyaml(self, ssheet):
        self.out_file = samplesheet.csv2yaml(ssheet)
        assert os.path.exists(self.out_file)
        with open(self.out_file) as in_handle:
            info = yaml.load(in_handle)
        return info
