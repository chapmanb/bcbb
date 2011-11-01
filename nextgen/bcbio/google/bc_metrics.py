#!/usr/bin/env python
"""Functions for getting barcode statistics from demultiplexing"""

import os
import re
import copy
from bcbio.utils import UnicodeReader
import bcbio.google.connection
import bcbio.google.document
import bcbio.google.spreadsheet
from bcbio.google import (_from_unicode,_to_unicode)
from bcbio.pipeline import log


# The structure of the demultiplex result
BARCODE_STATS_HEADER = [
                 ['Project name','project_name'],
                 ['Date','date'],
                 ['Flowcell','flowcell_id'],
                 ['Lane','lane'],
                 ['Description','description'],
                 ['Sample name','sample_name'],
                 ['Internal barcode index','barcode_id'],
                 ['Barcode name','name'],
                 ['Barcode sequence','sequence'],
                 ['Barcode type','barcode_type'],
                 ['Demultiplexed read count','barcode_read_count'],
                 ['Demultiplexed read count (millions)','barcode_read_count_millions'],
                 ['Comment','comment']
                ]
  
# The structure of the sequencing result
SEQUENCING_RESULT_HEADER = [
                 ['Sample name','sample_name'],
                 ['Run','run'],
                 ['Lane','lane'],
                 ['Sample count','read_count'],
                 ['Sample count (millions)','read_count_millions'],
                 ['Comment','comment'],
                 ['Pass','pass']
                ]
  
def _apply_filter(unfiltered,filter):
    """Remove rows whose column contents does not match the non-null contents in the filter columns"""
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
    
    # Get the GDocs demultiplex result file title
    gdocs_spreadsheet = gdocs.get("gdocs_dmplx_file",None)
    if not gdocs_spreadsheet:
        log.warn("Could not find Google Docs demultiplex results file title in config. No demultiplex counts were written to Google Docs")
        return
    
    # Get the account credentials
    encoded_credentials = ""
    encoded_credentials_file = gdocs.get("gdocs_credentials",None)
    if not encoded_credentials_file:
        log.warn("Could not find Google Docs account credentials. No demultiplex report was written")
        return
    # Check if the credentials file exists
    if not os.path.exists(encoded_credentials_file):
        log.warn("The Google Docs credentials file could not be found. No demultiplex data was written")
        return
    with open(encoded_credentials_file) as fh:
        encoded_credentials = fh.read().strip()
    
    # Get the barcode statistics. Get a deep copy of the run_info since we will modify it
    bc_metrics = get_bc_stats(fc_date,fc_name,work_dir,copy.deepcopy(run_info))
    
    # Upload the data
    write_run_report_to_gdocs(fc_date,fc_name,bc_metrics,gdocs_spreadsheet,encoded_credentials)
    
    # Get the projects parent folder
    projects_folder = gdocs.get("gdocs_projects_folder",None)
    
    # Write the bc project summary report
    if projects_folder:
        write_project_report_to_gdocs(fc_date,fc_name,bc_metrics,encoded_credentials,projects_folder)


def format_project_name(unformated_name):
    """Make the project name adhere to the formatting convention"""
    regexp = r'^(.+?)_(\d{2})_(\d{2})(.*)$'
    m = re.match(regexp,unformated_name)
    if not m or len(m.groups()) < 3:
        return unformated_name
    
    name = m.group(1).strip()
    year = m.group(2).strip()
    month = m.group(3).strip()
    suffix = m.group(4).strip()
   
    # Replace any non-period delimiters
    delimiter = "_"
    p = re.compile('(_)')
    name = p.sub('.',name)
    
    # Format the name
    project_name = "%s_%s_%s%s" % (name,year,month,suffix)
    return project_name
   
def get_bc_stats(fc_date, fc_name, work_dir, run_info):
    """Get a data structure with the run info coupled with the results from barcode demultiplexing"""
    bc_stats = []
    for lane_run_info in run_info.get("details",[]):
        lane_bc_stats = {}
        lane_id = str(lane_run_info['lane'])
        bc_dir = os.path.join(work_dir,"%s_%s_%s_barcode" % (lane_id,fc_date,fc_name))
        bc_file = os.path.join(bc_dir,"%s_%s_%s_bc.metrics" % (lane_id,fc_date,fc_name))
        if os.path.exists(bc_file):
            with open(bc_file) as bch:
                csvr = UnicodeReader(bch,dialect='excel-tab')
                for row in csvr:
                    lane_bc_stats[str(row[0])] = int(row[1])
        bc_stats.append(_merge_bc_stats(lane_run_info,lane_bc_stats,fc_date,fc_name))
    
    return bc_stats

def get_project_name(description):
    """Parse out the project name from the lane description"""
    m = re.match(r'(?:.*\s+)?(\S+)',description,re.I)
    if m and len(m.groups()) > 0:
        return format_project_name(m.group(1).strip())
    return "N/A"
       
def _get_unique_project_names(rows):
    """Get the unique project names in a set of rows"""
    names = {}
    for row in rows:
        names[row[0]] = 1
    return names.keys()

def get_sample_name(barcode_name):
    """Extract the sample name by stripping the barcode index part of the sample description""" 
    regexp = r'^(.+?)[\.\-_]?ind?(?:ex)?[\.\-_]?\d+$'
    m = re.search(regexp,(barcode_name or ""),re.I)
    if not m or len(m.groups()) == 0:
        return barcode_name
    return m.group(1)
    
def get_spreadsheet(ssheet_title,encoded_credentials):
    """Connect to Google docs and get a spreadsheet"""
    
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
        return (None,None)
    
    log.info("Found spreadsheet matching the supplied title: '%s'" % (ssheet.title.text))
    
    return (client,ssheet)

def group_bc_stats(bc_metrics):
    """Collapse the barcode statistics into a data structure where the number of reads per project and sample are aggregated"""
    
    projects = {}
    for lane in bc_metrics:
        lane_project_name = lane.get("project_name","N/A")
        lane_name = str(lane.get("lane","N/A"))
        
        for bc in lane.get("multiplex",[]):
            # If a project name is specified for the sample, use that instead of the lane project
            project_name = bc.get("project_name",lane_project_name)
            # Get the project data already stored on this project
            project = projects.get(project_name,None)
            if not project:
                project = {}
                project["project_name"] = project_name
                project["samples"] = {}
                projects[project_name] = project
            samples = project["samples"]
            
            sample_count = bc.get("barcode_read_count",0)
            sample_name = bc.get("sample_name","N/A")
            # Get the sample info already stored
            sample = samples.get(sample_name,None)
            if not sample:
                sample = {}
                sample["sample_name"] = sample_name
                sample["read_count"] = 0 
                sample["lane"] = []
                samples[sample_name] = sample
            sample["lane"].append(lane_name)
                
            # Convert the count to int
            try:
                sample_count = int(sample_count)
                # Add up the read counts
                sample["read_count"] = sample["read_count"] + sample_count
            except ValueError:
                sample["read_count_na"] = 1
    
    # Reformat the data structure
    data = projects.values()
    for project in data:
        samples = project["samples"].values()
        for sample in samples:
            lanes = dict.fromkeys(sample["lane"]).keys()
            sample["lane"] = ",".join(lanes)
            if "read_count_na" in sample:
                sample["comment"] = unicode("Read count not available for some lanes")
        project["samples"] = samples
        
    return data
       
def _merge_bc_stats(lane_run_info,lane_bc_stats,fc_date="N/A",fc_name="N/A"):
    """Join the barcode statistics with the run meta data"""
    
    lane_info = dict(lane_run_info)
    
    # Parse the project name
    lane_info['project_name'] = get_project_name(lane_info.get("description",""))
    lane_info['date'] = lane_info.get("date",fc_date)
    lane_info['flowcell_id'] = lane_info.get("flowcell_id",fc_name)
    
    # Add a multiplex section if none exists
    if not 'multiplex' in lane_info:
        lane_info['multiplex'] = []
         
    multiplex = lane_info['multiplex']
    for bc in multiplex:
        bc_index = str(bc['barcode_id'])
        bc['sample_name'] = get_sample_name(bc['name'])
        # set the project name based on the sample description or, if not present, the lane description
        bc['project_name'] = get_project_name(bc.get('description',lane_info.get("description","")))
        bc_count = lane_bc_stats.get(bc_index,None)
        if bc_count:
            bc['barcode_read_count'] = bc_count
            del lane_bc_stats[bc_index]
        else:
            bc['barcode_read_count'] = "N.D."
    
    # Add entries for barcodes not specified in the configuration file
    for bc_index, bc_count in lane_bc_stats.items():
        bc = {'barcode_id': bc_index, 'barcode_read_count': bc_count}
        # In case the barcode index is 'unmatched', use this as the sample name as well
        if bc_index == 'unmatched':
            bc['sample_name'] = 'Unmatched'
        multiplex.append(bc)
        
    return lane_info
 
def _structure_to_list(structure):
    """Flatten all entries in the metrics data structure into a list of entry rows"""
    
    metrics_list = []
    for lane in structure:
        row = [""]*len(BARCODE_STATS_HEADER)
        for i in range(0,5):
            row[i] = lane.get(BARCODE_STATS_HEADER[i][1],"")

        for m in lane.get("multiplex",[]):
            # use the sample-specific project names if present
            row[0] = m.get(BARCODE_STATS_HEADER[0][1],row[0])
            for i in range(5,len(BARCODE_STATS_HEADER)):
                row[i] = m.get(BARCODE_STATS_HEADER[i][1],"")
                
            metrics_list.append(list(row))
            
    return metrics_list
        
def write_project_report_to_gdocs(fc_date,fc_name,project_bc_metrics,encoded_credentials,gdocs_folder=""):
    """Upload the sample read distribution for a project to google docs"""
    
    # Create a client class which will make HTTP requests with Google Docs server.
    client = bcbio.google.spreadsheet.get_client(encoded_credentials)
    doc_client = bcbio.google.document.get_client(encoded_credentials)
    
    # Get a reference to the parent folder
    parent_folder = bcbio.google.document.get_folder(doc_client,gdocs_folder)
    
    # Group the barcode data by project
    grouped = group_bc_stats(project_bc_metrics)
    
    # Loop over the projects and write the project summary for each
    for pdata in grouped:
        
        project_name = pdata.get("project_name","")
        ssheet_title = project_name + "_sequencing_results"
        ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
        if not ssheet:
            bcbio.google.document.add_spreadsheet(doc_client,ssheet_title)
            ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
    
        _write_project_report_to_gdocs(client,ssheet,fc_date,fc_name,pdata)
        _write_project_report_summary_to_gdocs(client,ssheet)
        
        # Just to make it look a bit nicer, remove the default 'Sheet1' worksheet
        wsheet = bcbio.google.spreadsheet.get_worksheet(client,ssheet,'Sheet 1')
        if wsheet:
            client.DeleteWorksheet(wsheet)
            
        folder_name = project_name
        folder = bcbio.google.document.get_folder(doc_client,folder_name)
        if not folder:
            log.info("creating folder '%s'" % _from_unicode(folder_name))
            folder = bcbio.google.document.add_folder(doc_client,folder_name,parent_folder)
            
        ssheet = bcbio.google.document.move_to_folder(doc_client,ssheet,folder)
        log.info("'%s' spreadsheet written to folder '%s'" % (_from_unicode(ssheet.title.text),_from_unicode(folder_name)))
        

def _write_project_report_to_gdocs(client, ssheet, fc_date, fc_name, project_data):

    # Get the spreadsheet if it exists
    # Otherwise, create it
    wsheet_title = "%s_%s" % (fc_date,fc_name)
    
    # Flatten the project_data structure into a list
    rows = []
    for sample in project_data["samples"]:
        scount = int(sample["read_count"])
        mcount = round(scount/1000000.,2)
        row = (sample["sample_name"],"%s_%s" % (fc_date,fc_name),sample["lane"],unicode(scount),unicode(mcount),sample.get("comment",""),"")
        rows.append(row)
    
    # Write the data to the worksheet
    return _write_to_worksheet(client,ssheet,wsheet_title,rows,SEQUENCING_RESULT_HEADER,False)
    
def _write_project_report_summary_to_gdocs(client, ssheet):
    """Summarize the data from the worksheets and write them to a "Summary" worksheet"""
    
    # Summary data
    summary_data = {}
    
    # Get the list of worksheets in the spreadsheet
    wsheet_feed = bcbio.google.spreadsheet.get_worksheets_feed(client,ssheet)
    # Loop over the worksheets and parse the data from the ones that contain flowcell data
    for wsheet in wsheet_feed.entry:
        if not re.match(r'^\d{2}[01]\d[0-3]\d_\S+XX$',wsheet.title.text):
            continue
        wsheet_data = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet,'2')
        # Add the results from the worksheet to the summarized data    
        for wsheet_row in wsheet_data:
            sample_name = wsheet_row[0]
            summary_row = summary_data.get(sample_name,None)
            if not summary_row:
                summary_row = [""]*len(wsheet_row)
                summary_row[0] = sample_name
                for i in range(1,len(wsheet_row)):
                    summary_row[i] = []
                summary_data[sample_name] = summary_row 
            for i in range(1,len(wsheet_row)):
                summary_row[i].append(wsheet_row[i])
    
    # Concatenate the data in the summary structure
    for sample in summary_data.values():
        # Concatenate the non-number columns using ';' as separator
        for i in (1,2,5,6):
            sample[i] = ";".join(sample[i])
        
        # Sum up the read counts
        counts = sample[3]
        sum = 0
        for count in counts:
            sum += int(count)
        # Count the millions
        msum = round(sum/1000000.,2)
        sample[3] = unicode(sum)
        sample[4] = unicode(msum)
        
    # Write the summary data to the worksheet
    wsheet_title = "Summary"
    return _write_to_worksheet(client,ssheet,wsheet_title,summary_data.values(),SEQUENCING_RESULT_HEADER,False)
            

def write_run_report_to_gdocs(fc_date, fc_name, bc_metrics, ssheet_title, encoded_credentials, wsheet_title=None, append=False, split_project=False):
    """Upload the barcode read distribution for a run to google docs"""
    
    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title,encoded_credentials)
    if not client or not ssheet:
        return False
    
    # Convert the bc_metrics data structure into a flat list
    rows = _structure_to_list(bc_metrics)
    
    # Get the projects in the run
    projects = _get_unique_project_names(rows)
    log.info("The run contains data from: '%s'" % "', '".join(projects))
    
    # Calculate the number of million reads for convenience
    brci = -1
    brcmi = -1
    for i,head in enumerate(BARCODE_STATS_HEADER):
        if head[1] == 'barcode_read_count':
            brci = i
        elif head[1] == 'barcode_read_count_millions':
            brcmi = i
    if brci >= 0 and brcmi >= 0:
        for row in rows:
            try:
                row[brcmi] = unicode(round(int(row[brci])/1000000.,2))
            except ValueError:
                pass
    
    # If we will split the worksheet by project, use the project names as worksheet titles
    success = True
    if split_project:
        # Filter away the irrelevent project entries and write the remaining to the appropriate worksheet
        for wsheet_title in projects:
            success &= _write_to_worksheet(client,ssheet,wsheet_title,_apply_filter(rows,[wsheet_title]),BARCODE_STATS_HEADER,append)
            
    # Else, set the default title of the worksheet to be a string of concatenated date and flowcell id
    else:
        if wsheet_title is None:
            wsheet_title = "%s_%s" % (fc_date,fc_name)
        success &= _write_to_worksheet(client,ssheet,wsheet_title,rows,BARCODE_STATS_HEADER,append)

    return success

def _write_to_worksheet(client,ssheet,wsheet_title,rows,header,append):
    """Generic method to write a set of rows to a worksheet on google docs"""
    
    # Convert the worksheet title to unicode
    wsheet_title = _to_unicode(wsheet_title)
    
    # Add a new worksheet, possibly appending or replacing a pre-existing worksheet according to the append-flag
    wsheet = bcbio.google.spreadsheet.add_worksheet(client,ssheet,wsheet_title,len(rows)+1,len(header),append)
    if wsheet is None:
        log.info("Could not add a worksheet '%s' to spreadsheet '%s'" % (wsheet_title,ssheet.title.text))
        return False
    
    # Write the data to the worksheet
    log.info("Adding data to the '%s' worksheet" % (wsheet_title))
    return bcbio.google.spreadsheet.write_rows(client,ssheet,wsheet,[col_header[0] for col_header in header],rows)
     
    
    