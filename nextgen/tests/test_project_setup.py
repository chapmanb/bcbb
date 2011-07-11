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

class ProjectSetupTest(unittest.TestCase):
    """Test project setup"""

    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        self.fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        delivery_dir = os.path.join(self.proj_dir, "data", "110106_FC70BUKAAXX")
        if os.path.exists(delivery_dir):
            shutil.rmtree(delivery_dir)
        if os.path.exists(os.path.join(self.proj_dir, "intermediate")):
            shutil.rmtree(os.path.join(self.proj_dir, "intermediate"))
        self._deliver_data()

    def _deliver_data(self):
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.fcdir, self.proj_dir,
              "--flowcell_alias=20000101A_hiseq2000"]
        subprocess.check_call(cl)

    def test_project_setup(self):
        """Test project setup"""
        cl = ["setup_project_files.py",
              os.path.join(self.proj_dir, "data", "20000101A_hiseq2000", "project_run_info.yaml"),
              "20000101A_hiseq2000",
              "--project_dir=%s" %(self.proj_dir)]
        subprocess.check_call(cl)
