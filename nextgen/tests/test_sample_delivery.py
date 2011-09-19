"""
Test sample delivery
"""
import os
import sys
import subprocess
import unittest
import shutil
from string import Template

class SampleDeliveryTest(unittest.TestCase):
    """Deliver samples from bcbio-processed data to John Doe"""
    
    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        self.fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        delivery_dir = os.path.join(self.proj_dir, "data", "110106_FC70BUKAAXX")
        self._install_config_data()
        if os.path.exists(delivery_dir):
            shutil.rmtree(delivery_dir)
        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)

    def _install_config_data(self):
        loc_files = ['bowtie_indices.loc', 'bwa_index.loc', 'sam_fa_indices.loc']
        tooldir = os.path.join(self.file_dir, "config", "tool-data")

        d = {'genomedir':os.path.join(os.path.dirname(__file__), "data", "genomes")}
        if not os.path.exists(tooldir):
            os.makedirs(tooldir)
        for lf in loc_files:
            outfile = os.path.join(tooldir, lf)
            if not os.path.exists(outfile):
                with open(os.path.join(self.file_dir, "templates", "tool-data", lf)) as in_handle:
                    tmpl = Template(in_handle.read())
                with open(outfile, "w") as out_handle:
                    out_handle.write(tmpl.safe_substitute(d))
                          


    def test_deliver_data(self):
        """Test data delivery"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.fcdir, self.proj_dir]
        subprocess.check_call(cl)

    def test_dry_run(self):
        """Test dry run: don't do anything"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.fcdir, self.proj_dir, "-n"]
        subprocess.check_call(cl)

    def test_write_run_info_only(self):
        """Test writing of pruned run info file"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.fcdir, self.proj_dir, "-I"]
        subprocess.check_call(cl)
