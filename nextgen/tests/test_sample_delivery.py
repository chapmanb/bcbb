"""
Test sample delivery
"""
import os
import sys
import subprocess
import unittest
import shutil

try:
    import bcbio
except:
    raise ImportError("Module bcbio required to run sample based analysis. Make sure to run python setup.py develop so bcbio.__path__ is defined.")

class SampleDeliveryTest(unittest.TestCase):
    """Deliver samples from bcbio-processed data to John Doe"""
    
    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        self.bcbio_tests_dir = os.path.join(bcbio.__path__[0], os.pardir, "tests")
        self.bcbio_fcdir = os.path.join(self.bcbio_tests_dir, "test_automated_output")
        if os.path.exists(os.path.join(self.proj_dir, "data")):
            shutil.rmtree(os.path.join(self.proj_dir, "data"))

    def test_deliver_data(self):
        """Test data delivery"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.bcbio_fcdir, self.proj_dir]
        subprocess.check_call(cl)

    def test_dry_run(self):
        """Test dry run: don't do anything"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.bcbio_fcdir, self.proj_dir, "-n"]
        subprocess.check_call(cl)

    def test_write_run_info_only(self):
        """Test writing of pruned run info file"""
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.bcbio_fcdir, self.proj_dir, "-I"]
        subprocess.check_call(cl)
