"""
Test sample delivery
"""
import os
import sys
import subprocess
import unittest
import shutil
import yaml

class SampleDeliveryTest(unittest.TestCase):
    """Deliver samples from bcbio-processed data to John Doe"""
    
    def setUp(self):
        self.fcid = "110106_FC70BUKAAXX"
        self.file_dir = os.path.dirname(os.path.abspath(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects")
        self.fcdir = os.path.join(self.file_dir, "test_automated_output")
        self.ngsdata_dir  = os.path.join(self.file_dir, "ngsdata")
        self.archive_base_dir  = os.path.join(self.ngsdata_dir, "archive")
        self.analysis_base_dir  = os.path.join(self.ngsdata_dir, "analysis")
        self.post_process = os.path.join(self.analysis_base_dir, self.fcid, "post_process.yaml")
        if os.path.exists(self.proj_dir):
            shutil.rmtree(self.proj_dir)
        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)
        # if os.path.exists(self.ngsdata_dir):
        #     shutil.rmtree(self.ngsdata_dir)
        if not os.path.exists(self.analysis_base_dir):
            os.makedirs(self.analysis_base_dir)
            src = self.fcdir
            dest = os.path.abspath(os.path.join(self.analysis_base_dir, self.fcid))
            shutil.copytree(src, dest)
        if not os.path.exists(self.archive_base_dir):
            os.makedirs(os.path.join(self.archive_base_dir, self.fcid))
            src = os.path.abspath(os.path.join(self.file_dir, "data", "automated", "run_info-project.yaml"))
            dest = os.path.join(self.archive_base_dir, self.fcid, "run_info.yaml")
            os.symlink(src, dest)
        with open(os.path.join(self.file_dir, "data", "automated", "post_process.yaml"), "r") as fh:
            post_process = yaml.load(fh)
        post_process["analysis"]["store_dir"] = os.path.join(self.archive_base_dir)
        post_process["analysis"]["base_dir"] = os.path.join(self.analysis_base_dir)
        with open(os.path.join(self.analysis_base_dir, self.fcid, "post_process.yaml"), "w") as fh:
            yaml.dump(post_process, stream=fh)

    def test_a_dry_run(self):
        """Test dry run: don't do anything"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01",
              "--dry_run"]
        subprocess.check_call(cl)

    def test_a_write_run_info_only(self):
	"""Test writing of pruned run info file"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01",
              "--only_install_run_info"]
        subprocess.check_call(cl)

    def test_b_deliver_data_copy(self):
        """Test data delivery by default method copy"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01"]
        subprocess.check_call(cl)

    def test_deliver_data(self):
        """Test data delivery by moving"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir, 
              "--project_desc=%s" % "J.Doe_00_01", "--move_data"]
        subprocess.check_call(cl)

