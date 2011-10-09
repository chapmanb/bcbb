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
        self.file_dir = os.path.dirname(os.path.abspath(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects")
        self.fcdir = os.path.join(self.file_dir, "test_automated_output")
        self.archive_base_dir  = os.path.join(self.file_dir, "scilife", "archive")
        self.analysis_base_dir  = os.path.join(self.file_dir, "scilife", "analysis")
        self.scilife_dir  = os.path.join(self.file_dir, "scilife")

        self._install_config_data()
        if os.path.exists(self.proj_dir):
            shutil.rmtree(self.proj_dir)
        os.makedirs(self.proj_dir)
        if os.path.exists(self.scilife_dir):
            shutil.rmtree(self.scilife_dir)
        if not os.path.exists(self.analysis_base_dir):
            os.makedirs(self.analysis_base_dir)
            src = self.fcdir
            dest = os.path.abspath(os.path.join(self.analysis_base_dir, "110106_FC70BUKAAXX"))
            shutil.copytree(src, dest)
        if not os.path.exists(self.archive_base_dir):
            os.makedirs(os.path.join(self.archive_base_dir, "110106_FC70BUKAAXX"))
            src = os.path.abspath(os.path.join(self.file_dir, "data", "automated", "run_info.yaml"))
            src = os.path.abspath(os.path.join(self.file_dir, "templates", "run_info.yaml"))
            dest = os.path.join(self.archive_base_dir, "110106_FC70BUKAAXX", "run_info.yaml")
            os.symlink(src, dest)


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
              "110106_FC70BUKAAXX", "j_doe_00_01",
              "--analysis_base_dir=%s" % self.analysis_base_dir,
              "--archive_base_dir=%s" % self.archive_base_dir,
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01"]
        subprocess.check_call(cl)

    def test_dry_run(self):
        """Test dry run: don't do anything"""
        self.setUp()
        cl = ["sample_delivery.py",
              "110106_FC70BUKAAXX", "j_doe_00_01",
              "--analysis_base_dir=%s" % self.analysis_base_dir,
              "--archive_base_dir=%s" % self.archive_base_dir,
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01",
              "--dry_run"]
        subprocess.check_call(cl)

    def test_write_run_info_only(self):
	"""Test writing of pruned run info file"""
        self.setUp()
        cl = ["sample_delivery.py",
              "110106_FC70BUKAAXX", "j_doe_00_01",
              "--analysis_base_dir=%s" % self.analysis_base_dir,
              "--archive_base_dir=%s" % self.archive_base_dir,
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01",
              "--only_install_run_info"]
        subprocess.check_call(cl)
