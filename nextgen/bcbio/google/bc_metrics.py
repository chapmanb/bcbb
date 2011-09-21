#!/usr/bin/env python
"""Functions for getting barcode statistics from demultiplexing"""

import os
import re

from bcbio.utils import UnicodeReader
from bcbio.pipeline import log
import bcbio.google.connection
import bcbio.google.spreadsheet

# The structure of the demultiplex result file on the form {title: column index}
COLUMN_HEADER = [
                 ['Project name','project_name'],
                 ['Date','date'],
                 ['Flowcell','flowcell_id'],
                 ['Lane','lane'],
                 ['Description','description'],
                 ['Internal barcode index','barcode_id'],
                 ['Barcode name','name'],
                 ['Barcode sequence','sequence'],
                 ['Barcode type','barcode_type'],
                 ['Demultiplexed read count','barcode_read_count'],
                 ['Comment','comment']
                ]
  
def _apply_filter(unfiltered,filter):
    
    filtered = []
    for entry in unfiltered:
        passed = True
        for i,f in enumerate(filter):
            if f and entry[i] != f:
                passed = False
                break
        if passed:
            filtered.append(entry)
    
    return filtered

def create_bc_report_on_gdocs(fc_date, fc_name, work_dir, run_info, config):
    """Get the barcode read distribution for a run and upload to google docs"""
    
    # Get the required parameters from the post_process.yaml configuration file
    gdocs = config.get("gdocs_upload",None)
    if not gdocs:
        log.info("No GDocs upload section specified in config file, will not upload demultiplex data")
        return
    
    encoded_credentials = gdocs.get("gdocs_credentials",None)
    
    # Get the GDocs demultiplex result file title
    gdocs_spreadsheet = gdocs.get("gdocs_dmplx_file",None)
    if not gdocs_spreadsheet:
        log.warn("Could not find Google Docs demultiplex results file title in config. No demultiplex counts were written to Google Docs")
        return
    
    # Get the account credentials
    encoded_credentials = gdocs.get("gdocs_credentials",None)
    if not encoded_credentials:
        log.warn("Could not find Google Docs account credentials. No demultiplex report was written")
        return
    
    # Get the barcode statistics
    bc_metrics = get_bc_stats(fc_date,fc_name,work_dir,run_info)
    
    # Upload the data
    write_bc_report_to_gdocs(fc_date,fc_name,bc_metrics,gdocs_spreadsheet,encoded_credentials)
   
def _from_unicode(unistr,encoding='utf-8'):
    if isinstance(unistr,unicode):
        unistr = unistr.encode(encoding)
    return unistr

   
def get_bc_stats(fc_date, fc_name, work_dir, run_info):
    """Get a data structure with the run info coupled with the results from barcode demultiplexing"""
    
    bc_stats = []
    for lane_run_info in run_info:
        lane_bc_stats = {}
        lane_id = unicode(lane_run_info['lane'])
        bc_dir = os.path.join(work_dir,"%s_%s_%s_barcode" % (lane_id,fc_date,fc_name))
        bc_file = os.path.join(bc_dir,"%s_%s_%s_bc.metrics" % (lane_id,fc_date,fc_name))
        if os.path.exists(bc_file):
            with open(bc_file) as bch:
                csvr = UnicodeReader(bch,dialect='excel-tab')
                for row in csvr:
                    lane_bc_stats[row[0]] = int(row[1])
        bc_stats.append(merge_bc_stats(lane_run_info,lane_bc_stats,fc_date,fc_name))
        
    return bc_stats

def get_project_name(description):
    
    m = re.match(r'Lane\s+\S+\s*,\s*(.*)',description,re.I)
    if m and len(m.groups()) > 0:
        return m.group(1).strip()
    return "N/A"
       
def _get_project_names(rows):
    names = {}
    for row in rows:
        names[row[0]] = 1
    return names.keys()
                          
def merge_bc_stats(lane_run_info,lane_bc_stats,fc_date="N/A",fc_name="N/A"):
    """Get a data structure with lane and sample meta data merged with barcode statistics
    
    --lane_run_info - the data structure for a lane from a parsed run_info.yaml configuration file
    --lane_bc_stats - a dictionary with the internal barcode index as keys and the corresponding read counts as values"""
    
    lane_info = dict(lane_run_info)
    
    # Parse the project name based on the assumption that the description is on the form "[Lane id], [Project name]"
    lane_info['project_name'] = get_project_name(lane_info.get("description",""))
    lane_info['date'] = lane_info.get("date",fc_date)
    lane_info['flowcell_id'] = lane_info.get("flowcell_id",fc_name)
    
    # Add a multiplex section if none exists
    if not 'multiplex' in lane_info:
        lane_info['multiplex'] = []
         
    multiplex = lane_info['multiplex']
    for bc in multiplex:
        bc_index = unicode(bc['barcode_id'])
        bc_count = lane_bc_stats.get(bc_index,None)
        if bc_count:
            bc['barcode_read_count'] = bc_count
            del lane_bc_stats[bc_index]
        else:
            bc['barcode_read_count'] = "N.D."
    
    # Add entries for barcodes not specified in the configuration file
    for bc_index, bc_count in lane_bc_stats.items():
        multiplex.append({'barcode_id': bc_index, 'barcode_read_count': bc_count})
        
    return lane_info
 
def _structure_to_list(structure):
    """Flatten all entries in the metrics data structure into a list of entry rows"""
    
    metrics_list = []
    for lane in structure:
        row = [""]*len(COLUMN_HEADER)
        for i in range(0,5):
            row[i] = lane.get(COLUMN_HEADER[i][1],"")

        for m in lane.get("multiplex",[]):
            for i in range(5,len(COLUMN_HEADER)):
                row[i] = m.get(COLUMN_HEADER[i][1],"")
                
            metrics_list.append(list(row))
            
    return metrics_list

def _to_unicode(str,encoding='utf-8'):
    if isinstance(str,basestring):
        if not isinstance(str,unicode):
            str = unicode(str,encoding)
    return str

def write_bc_report_to_gdocs(fc_date, fc_name, bc_metrics, ssheet_title, encoded_credentials, wsheet_title=None, append=False, split_project=False):
    """Upload the barcode read distribution to google docs"""
    
    # Convert the spreadsheet title to unicode
    ssheet_title = _to_unicode(ssheet_title)
    
    # Create a client class which will make HTTP requests with Google Docs server.
    client = bcbio.google.spreadsheet.get_client()
    bcbio.google.connection.authenticate(client,encoded_credentials)
    
    # Locate the spreadsheet
    ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
    
    # Check that we got a result back
    if not ssheet:
        log.warn("No document with specified title '%s' found in GoogleDocs repository" % ssheet_title)
        return False
    
    log.info("Found spreadsheet matching the supplied title: '%s'" % (ssheet.title.text))
    
    # Convert the bc_metrics data structure into a flat list
    rows = _structure_to_list(bc_metrics)
    
    # Get the projects in the run
    projects = _get_project_names(rows)
    log.info("The run contains data from: '%s'" % "', '".join(projects))
    
    # If we will split the worksheet by project, use the project names as worksheet titles
    success = True
    if split_project:
        # Filter away the irrelevent project entries and write the remaining to the appropriate worksheet
        for wsheet_title in projects:
            success &= _write_to_worksheet(client,ssheet,wsheet_title,_apply_filter(rows,[wsheet_title]),append)
            
    # Else, set the default title of the worksheet to be a string of concatenated date and flowcell id
    else:
        if wsheet_title is None:
            wsheet_title = "%s_%s" % (fc_date,fc_name)
        success &= _write_to_worksheet(client,ssheet,wsheet_title,rows,append)

    return success

def _write_to_worksheet(client,ssheet,wsheet_title,rows,append):
    
    # Convert the worksheet title to unicode
    wsheet_title = _to_unicode(wsheet_title)
    
    # Add a new worksheet, possibly appending or replacing a pre-existing worksheet according to the append-flag
    wsheet = bcbio.google.spreadsheet.add_worksheet(client,ssheet,wsheet_title,len(rows)+1,len(COLUMN_HEADER),append)
    if wsheet is None:
        log.info("Could not add a worksheet '%s' to spreadsheet '%s'" % (wsheet_title,ssheet.title.text))
        return False
    
    # Write the data to the worksheet
    log.info("Adding data to the '%s' worksheet" % (wsheet_title))
    return bcbio.google.spreadsheet.write_rows(client,ssheet,wsheet,[col_header[0] for col_header in COLUMN_HEADER],rows)
 
    
    