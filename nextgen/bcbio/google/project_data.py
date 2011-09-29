#!/usr/bin/env python
"""Functions for doing google docs-related stuff centered around a project"""

import bcbio.google.connection
import bcbio.google.spreadsheet
from bcbio.pipeline import log


def get_uppnex_project_id(project_name,config):
    """Attempt to fetch the Uppnex id associated with a project"""
    
    PROJECT_NAME_COLUMN = 'Project name'
    UPPNEXID_COLUMN = 'Uppnex ID'
    
    # Get the name of the spreadsheet where uppnex ids can be found
    gdocs_config = config.get("gdocs_upload",{})
    ssheet_title = gdocs_config.get("projects_spreadsheet",None)
    wsheet_title = gdocs_config.get("projects_worksheet",None)
    if not ssheet_title or not wsheet_title:
        log.warn("The names of the projects spreadsheet and worksheet on Google Docs could not be found. The Uppnex ID could not be determined.")
        return 'N/A'
    
    # Get the account credentials
    encoded_credentials = gdocs_config.get("gdocs_credentials",None)
    if not encoded_credentials:
        log.warn("Could not find Google Docs account credentials. The Uppnex ID could not be determined.")
        return 'N/A'
    
    # Connect to the spread- and worksheet
    client = bcbio.google.spreadsheet.get_client()
    bcbio.google.connection.authenticate(client,encoded_credentials)
    ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
    if not ssheet:
        log.warn("Could not connect to %s on Google Docs. The Uppnex ID could not be determined." % ssheet_title)
        return 'N/A'
    wsheet = bcbio.google.spreadsheet.get_worksheet(client,ssheet,wsheet_title)
    if not wsheet:
        log.warn("Could not locate %s in %s. The Uppnex ID could not be determined." % (wsheet_title,ssheet_title))
        return 'N/A'
    
    # Get the Uppnex id for the project
    uppnex_id_column = bcbio.google.spreadsheet.get_column(client,ssheet,wsheet,UPPNEXID_COLUMN,{PROJECT_NAME_COLUMN: project_name})
    if len(uppnex_id_column) == 0:
        return 'N/A'
    
    return ", ".join(uppnex_id_column)
   