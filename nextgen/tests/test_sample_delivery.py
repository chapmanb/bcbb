"""
Test sample delivery
"""
import os
import sys
import subprocess
import unittest
import shutil
import yaml
import contextlib
from bcbio.pipeline.config_loader import load_config

@contextlib.contextmanager
def make_workdir(link=True):
    dirname = os.path.join(os.path.dirname(__file__), "110106_FC70BUKAAXX")
    if os.path.exists(dirname):
        if os.path.islink(dirname):
            os.remove(dirname)
        else:
            shutil.rmtree(dirname)
    src = os.path.join(os.path.dirname(__file__), "test_automated_output")
    orig_dir = os.getcwd()
    if link:
        os.symlink(src, dirname)
    else:
        shutil.copytree(src, dirname)
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)
    
class SampleDeliveryTest(unittest.TestCase):
    """Deliver samples from bcbio-processed data to John Doe"""

    def setUp(self):
        self.fcid = "110106_FC70BUKAAXX"
        self.file_dir = os.path.dirname(os.path.abspath(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects")
        self.fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        self.archive_base_dir  = os.path.join(self.file_dir)
        self.analysis_base_dir  = os.path.join(self.file_dir)

        self.post_process = os.path.join(self.analysis_base_dir, self.fcid, "post_process.yaml")
        if os.path.exists(self.proj_dir):
            shutil.rmtree(self.proj_dir)
        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "run_info.yaml")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "run_info-project.yaml"), os.path.join(self.file_dir, "test_automated_output", "run_info.yaml"))
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "tool-data")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "tool-data"), os.path.join(self.file_dir, "test_automated_output", "tool-data"))

        post_process = load_config(os.path.join(self.file_dir, "data", "automated", "post_process.yaml"))
        post_process["analysis"]["store_dir"] = os.path.join(self.archive_base_dir)
        post_process["analysis"]["base_dir"] = os.path.join(self.analysis_base_dir)
        with open(os.path.join(self.file_dir, "test_automated_output", "post_process.yaml"), "w") as fh:
            yaml.dump(post_process, stream=fh)

    def test_1_dry_run(self):
        """Test dry run: don't do anything"""
        self.setUp()
        with make_workdir():
            cl = ["sample_delivery.py",
                  self.post_process,
                  self.fcid, "j_doe_00_01",
                  "--project_base_dir=%s" % self.proj_dir,
                  "--project_desc=%s" % "J.Doe_00_01",
                  "--dry_run"]
            subprocess.check_call(cl)

    def test_2_write_run_info_only(self):
        """Test writing of pruned run info file"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir,
              "--project_desc=%s" % "J.Doe_00_01",
              "--only_install_run_info"]
        subprocess.check_call(cl)

    def test_3_deliver_data_copy(self):
        """Test data delivery by default method copy"""
        cl = ["sample_delivery.py",
              self.post_process,
              self.fcid, "j_doe_00_01",
              "--project_base_dir=%s" % self.proj_dir,
              "--project_desc=%s" % "J.Doe_00_01"]
        subprocess.check_call(cl)

    def test_4_deliver_data(self):
        """Test data delivery by moving"""
        with make_workdir(link=False):
            cl = ["sample_delivery.py",
                  self.post_process,
                  self.fcid, "j_doe_00_01",
                  "--project_base_dir=%s" % self.proj_dir, 
                  "--project_desc=%s" % "J.Doe_00_01", "--move_data"]
            subprocess.check_call(cl)

