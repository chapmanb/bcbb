"""This directory is setup with configurations to run the main functional test.

It exercises a samplebased analysis pipeline on a smaller subset of data, as implemented at SciLife.
"""
import os
import sys
import subprocess
import unittest
import shutil
import contextlib
import glob
import yaml

@contextlib.contextmanager
def workdir():
    dirname = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "intermediate", "nobackup", "20000101A_hiseq2000")
    # if os.path.exists(dirname):
    #     shutil.rmtree(dirname)
    # os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)

class SampleBasedAnalysisTest(unittest.TestCase):
    """Setup a sample based scilife analysis
    """
    def setUp(self):
        self.fcid = "110106_FC70BUKAAXX"
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects")
        self.fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        self.archive_base_dir  = os.path.join(self.file_dir)
        self.analysis_base_dir = os.path.join(self.file_dir)

        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)
        if not os.path.exists(os.path.join(self.file_dir, self.fcid)):
            os.symlink(os.path.join(self.file_dir, "test_automated_output"), os.path.join(self.file_dir, self.fcid))
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "run_info.yaml")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "run_info-project.yaml"), os.path.join(self.file_dir, "test_automated_output", "run_info.yaml"))
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "tool-data")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "tool-data"), os.path.join(self.file_dir, "test_automated_output", "tool-data"))

        # Post_process.yaml
        with open(os.path.join(self.file_dir, "data", "automated", "post_process.yaml"), "r") as fh:
            post_process = yaml.load(fh)
        post_process["analysis"]["store_dir"] = os.path.join(self.archive_base_dir)
        post_process["analysis"]["base_dir"] = os.path.join(self.analysis_base_dir)
        post_process["algorithm"]["snpcall"] = "true"
        post_process["algorithm"]["dbsnp"] = os.path.join("data", "genomes", "hg19", "variation", "dbsnp_132.vcf")
        with open(os.path.join(self.analysis_base_dir, self.fcid, "post_process.yaml"), "w") as fh:
            yaml.dump(post_process, stream=fh)
        self._deliver_data()
        
    def _deliver_data(self):
        print "Delivering data"
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, self.fcid, "post_process.yaml"),
              self.fcid, "j_doe_00_01",
              "--project_desc=%s" % "J.Doe_00_01",
              "--project_base_dir=%s" % self.proj_dir, 
              "--flowcell_alias=20000101A_hiseq2000"]
        subprocess.check_call(cl)
        print "Finished delivering data..."

    def test_run_samplebased_pipeline(self):
        """Test a sample based pipeline"""
        with workdir():
            cl = ["exome_pipeline.py",
                  os.path.join(self.analysis_base_dir, self.fcid, "post_process.yaml"),
                  os.path.join(self.proj_dir, "j_doe_00_01", "intermediate", "nobackup", "110106_FC70BUKAAXX"),
                  os.path.join(self.proj_dir, "j_doe_00_01", "data", "nobackup", "20000101A_hiseq2000", "project_run_info.yaml"),
                  "--project_dir=%s" %(self.proj_dir)]
            subprocess.check_call(cl)
