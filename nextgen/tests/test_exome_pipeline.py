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
from string import Template

try:
    import bcbio
except:
    raise ImportError("Module bcbio required to run sample based analysis. Make sure to run python setup.py develop so bcbio.__path__ is defined.")

@contextlib.contextmanager
def make_workdir():
    dirname = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "intermediate")
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
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
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        self.bcbio_tests_dir = os.path.join(bcbio.__path__[0], os.pardir, "tests")
        self.bcbio_fcdir = os.path.join(self.bcbio_tests_dir, "test_automated_output")
        self._install_project_config_files()
        if os.path.exists(os.path.join(self.proj_dir, "data")):
            shutil.rmtree(os.path.join(self.proj_dir, "data"))
        if os.path.exists(os.path.join(self.proj_dir, "intermediate")):
            shutil.rmtree(os.path.join(self.proj_dir, "intermediate"))
        self._deliver_data()
        self._setup_project()

    def _deliver_data(self):
        print "Delivering data"
        cl = ["sample_delivery.py",
              os.path.join(self.file_dir, "templates", "run_info.yaml"),
              "J.Doe_00_01", self.bcbio_fcdir, self.proj_dir,
              "--flowcell_alias=20000101A_hiseq2000"]
        subprocess.check_call(cl)
        print "Finished delivering data..."

    def _setup_project(self):
        print "setting up project"
        cl = ["setup_project_files.py",
              os.path.join(self.proj_dir, "data", "20000101A_hiseq2000", "project_run_info.yaml"),
              "20000101A_hiseq2000",
              "--project_dir=%s" %(self.proj_dir)]
        subprocess.check_call(cl)

        
    def _install_config_data(self):
        loc_files = ['bowtie_indices.loc', 'bwa_index.loc', 'sam_fa_indices.loc']
        tooldir = os.path.join(self.file_dir, "config", "tool-data")

        d = {'genomedir':os.path.join(self.bcbio_tests_dir, "data", "genomes")}
        if not os.path.exists(tooldir):
            os.makedirs(tooldir)
        for lf in loc_files:
            outfile = os.path.join(tooldir, lf)
            if not os.path.exists(outfile):
                with open(os.path.join(self.file_dir, "templates", "tool-data", lf)) as in_handle:
                    tmpl = Template(in_handle.read())
                with open(outfile, "w") as out_handle:
                    out_handle.write(tmpl.safe_substitute(d))
                          

    def _install_project_config_files(self): 
        """Install the project config files, inserting the correct paths to galaxy configuration files"""
        proj_conf = os.path.join(self.file_dir, "templates", "proj_conf.yaml")
        d = {'galaxy_config' : os.path.join(self.file_dir, "config", "universe_wsgi.ini"),
             'log_dir' : os.path.join(self.file_dir, "log")}
        with open(proj_conf) as in_handle:
            tmpl = Template(in_handle.read())
        print "Installing project configuration in " + self.proj_dir
        with open (os.path.join(self.proj_dir, "proj_conf.yaml"), "w") as out_handle:
            out_handle.write(tmpl.safe_substitute(d))

    def test_run_samplebased_pipeline(self):
        """Test a sample based pipeline"""
        cl = ["exome_pipeline.py",
              os.path.join(self.proj_dir, "proj_conf.yaml"),
              os.path.join(self.proj_dir, "data", "20000101A_hiseq2000"),
              os.path.join(self.proj_dir, "data", "20000101A_hiseq2000", "project_run_info.yaml"),
              "--project_dir=%s" %(self.proj_dir)]
        subprocess.check_call(cl)
