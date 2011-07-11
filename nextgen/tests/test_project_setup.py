"""
Test project setup

Once project fastq files have been delivered with sample_delivery.py they have 
to be relinked to comply with bcbio.
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

class ProjectSetupTest(unittest.TestCase):
    """Test project setup"""

    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        self.bcbio_tests_dir = os.path.join(bcbio.__path__[0], os.pardir, "tests")
        self.bcbio_fcdir = os.path.join(self.bcbio_tests_dir, "test_automated_output")
        if os.path.exists(os.path.join(self.proj_dir, "data")):
            shutil.rmtree(os.path.join(self.proj_dir, "data"))
        if os.path.exists(os.path.join(self.proj_dir, "intermediate")):
            shutil.rmtree(os.path.join(self.proj_dir, "intermediate"))
        self._deliver_data()

    def _deliver_data(self):
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.bcbio_fcdir, self.proj_dir,
              "--flowcell_alias=20000101A_hiseq2000"]
        subprocess.check_call(cl)

    def test_project_setup(self):
        """Test project setup"""
        cl = ["setup_project_files.py",
              os.path.join(self.proj_dir, "data", "20000101A_hiseq2000", "project_run_info.yaml"),
              "20000101A_hiseq2000",
              "--project_dir=%s" %(self.proj_dir)]
        subprocess.check_call(cl)
