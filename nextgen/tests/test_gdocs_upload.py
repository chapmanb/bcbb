"""A test that will attempt to generate demultiplex reports and upload to Google Docs
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import yaml
import random
from test_automated_analysis import make_workdir
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.utils import UnicodeWriter
from bcbio.google.bc_metrics import create_bc_report_on_gdocs

class GDocsUploadTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        make_workdir()
        # Make up some barcode numbers
        self.workdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")
        
        if os.path.exists(self.workdir):
            shutil.rmtree(self.workdir)
        os.makedirs(self.workdir)
        
        # Make up a bogus run name
        self.runname = "111014_SN0000_0000_AB0AAAACXX"
        
        # Parse the run_info
        run_info_file = os.path.join(self.data_dir, "run_info-gdocs.yaml")
        with open(run_info_file) as fh:
            self.run_info = yaml.load(fh)
            
        self._make_bc_metrics()
        

    def _make_bc_metrics(self):
        """Parses the run_info and generates lane folders and barcode metrics corresponding to the lanes and barcodes used"""
        fc_name, fc_date = get_flowcell_info(self.runname)
        barcode_dir_suffix = "_%s_%s_barcode" % (fc_date,fc_name)
        
        for lane in self.run_info:
            lane_name = str(lane['lane'])
            bc_dir = os.path.join(self.workdir,"%s%s" % (lane_name,barcode_dir_suffix))
            
            # Create the directory if it doesn't exist
            if not os.path.exists(bc_dir):      
                os.makedirs(bc_dir)
            
            # Create, or if it exists, append to the bc_metrics file
            bc_file = os.path.join(bc_dir,"%s_%s_%s_bc.metrics" % (lane_name,fc_date,fc_name))
            with open(bc_file,"a") as fh:
                bcw = UnicodeWriter(fh,dialect='excel-tab')
                
                # Loop over the barcodes and generate random read counts
                bcs = lane.get("multiplex",[])
                for bc in bcs:
                    bc_id = str(bc['barcode_id'])
                    bc_count = random.randint(1,10000000)
                    bcw.writerow([bc_id,bc_count])
                # Lastly write some unmatched counts, or in case no multiplex data was given, a 'trim' entry
                if len(bcs):
                    bcw.writerow(['unmatched',random.randint(1,10000000)])
                else:
                    bcw.writerow(['trim',random.randint(1,100000000)])
        
         
    def test_create_bc_report(self):
        """Create a demultiplex report and upload it to gdocs
        """
        fc_name, fc_date = get_flowcell_info(self.runname)
        
        # Parse the config
        config_file = os.path.join(self.data_dir, "post_process.yaml")
        with open(config_file) as fh:
            self.config = yaml.load(fh)
            
        create_bc_report_on_gdocs(fc_date, fc_name, self.workdir, {'details': self.run_info}, self.config)
        
