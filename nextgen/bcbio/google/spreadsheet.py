#!/usr/bin/env python
"""Wrapper functions around the python gdata spreadsheet api"""

import gdata.spreadsheet.service

def add_worksheet(client,ssheet,title,rows=0,cols=0,append=False):
    """Add a new worksheet with the specified title to the specified spreadsheet. 
    Will overwrite an existing worksheet with the same title unless append is True
    """
    
    # Check if a worksheet with the same title exists
    ws = get_worksheet(client,ssheet,title)
    if ws:
        # If we're appending, just return the first object in the feed
        if append:
            return ws
        
        # Otherwise, drop the existing worksheet
        client.DeleteWorksheet(ws)
    
    # Add the desired worksheet
    return client.AddWorksheet(title,rows,cols,get_key(ssheet))
    
def get_client():
    """Get a spreadsheet client"""
    # Create a client class which will make HTTP requests with Google Docs server.
    return gdata.spreadsheet.service.SpreadsheetsService()

def get_key(object):
    """Get the unique key identifier for the supplied object"""
    return object.id.text.split('/')[-1]

def _get_query(title,exact_match):
    """Get a query object for the supplied parameters"""
    
    p = {}
    if title:
        p['title'] = title
        if exact_match:
            p['title-exact'] = 'true'
        else:
            p['title-exact'] = 'false'
    
    return gdata.spreadsheet.service.DocumentQuery(params=p)

def get_spreadsheet(client,title):
    """Get an exact match for a spreadsheet"""
    feed = get_spreadsheets_feed(client,title,True)
    if len(feed.entry) == 0:
        return None
    return feed.entry[0]

def get_spreadsheets_feed(client, title=None, exact_match=False):
    """Get a feed of all available spreadsheets, optionally restricted by title"""
    
    # Create a query that restricts the spreadsheet feed to documents having the supplied title
    q = _get_query(title,exact_match)
    # Query the server for an Atom feed containing a list of your documents.
    return client.GetSpreadsheetsFeed(query=q)

def get_worksheet(client,ssheet,title):
    """Get an exact match for a worksheet within a spreadsheet"""
    feed = get_worksheets_feed(client,ssheet,title,True)
    if len(feed.entry) == 0:
        return None
    return feed.entry[0]

def get_worksheets_feed(client, ssheet, title=None, exact_match=False):
    """Get a feed of all worksheets in the supplied spreadsheet, optionally restricted by title"""
    
    # Create a query that restricts the spreadsheet feed to documents having the supplied title
    q = _get_query(title,exact_match)
    # Get the key for the spreadsheet
    k = get_key(ssheet)
    # Query the server for an Atom feed containing a list of your documents.
    return client.GetWorksheetsFeed(key=k,query=q)

def write_rows(client,ssheet,wsheet,header,rows):
    """Write the supplied data rows to the worksheet, using the supplied column headers"""
    
    # Get the keys
    ss_key = get_key(ssheet)
    ws_key = get_key(wsheet)
    
    try:
        # As a workaround for the InsertRow bugs with column names, just use single lowercase letters as column headers to start with
        for i in range(0,len(header)):
            client.UpdateCell(1,i+1,chr(97+i),ss_key,ws_key)
          
        # Iterate over the rows and add the data to the worksheet
        for row in rows:
            row_data = {}
            for i,value in enumerate(row):
                row_data[chr(97+i)] = unicode(value)
            client.InsertRow(row_data,ss_key,ws_key)

        # Lastly, substitute the one-letter header for the real string
        for i in range(0,len(header)):
            client.UpdateCell(1,i+1,header[i],ss_key,ws_key)
    except:
        return False
    
    return True
        

    