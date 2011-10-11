#!/usr/bin/env python
"""Wrapper functions around the python gdata spreadsheet api"""

import gdata.spreadsheet.service
import gdata.docs.service
from bcbio.google.connection import authenticate
from bcbio.google import (_from_unicode,_to_unicode)

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
    return client.AddWorksheet(_to_unicode(title),rows,cols,get_key(ssheet))

def column_count(wsheet):
    """Get the number of columns in the worksheet"""
    return int(wsheet.col_count.text)

def get_cell_content(client, ssheet, wsheet, row_start=0, col_start=0, row_end=0, col_end=0):
    """Get the text contents of the cells from the supplied spreadsheet and worksheet and from the specified cell range as a two-dimensional list"""
    
    if str(row_start) == '0':
        row_start = '1'
    if str(col_start) == '0':
        col_start = '1'
    if str(row_end) == '0':
        row_end = str(row_count(wsheet))
    if str(col_end) == '0':
        col_end = str(column_count(wsheet))
    
    feed = (get_cell_feed(client,ssheet,wsheet,row_start,col_start,row_end,col_end) or [])
    
    # Get the dimensions of the 2D-list
    rows = int(row_end) - int(row_start) + 1
    cols = int(col_end) - int(col_start) + 1
    content = []
    for i, cell in enumerate(feed.entry):
        r = i//cols
        c = i - r*cols
        if c == 0:
            row = []
            content.append(row)
        row.append(_to_unicode((cell.content.text or "")))
    
    return content
        

def get_cell_feed(client, ssheet, wsheet, row_start=0, col_start=0, row_end=0, col_end=0):
    """Get a cell feed from the supplied spreadsheet and worksheet and from the specified cell range"""
    
    if str(row_start) == '0':
        row_start = '1'
    if str(col_start) == '0':
        col_start = '1'
    if str(row_end) == '0':
        row_end = str(row_count(wsheet))
    if str(col_end) == '0':
        col_end = str(column_count(wsheet))
    
    p = {'min-row': str(row_start),
         'min-col': str(col_start),
         'max-row': str(row_end),
         'max-col': str(col_end),
         'return-empty': 'True'
         }
    query = gdata.spreadsheet.service.CellQuery(params=p)
    return client.GetCellsFeed(get_key(ssheet),get_key(wsheet),query=query)
    
def get_client(encoded_credentials=None):
    """Get a SpreadsheetsService client"""
    # Create a client class which will make HTTP requests with Google Docs server.
    client = gdata.spreadsheet.service.SpreadsheetsService()
    # If credentials were supplied, authenticate the client as well
    if encoded_credentials:
        authenticate(client,encoded_credentials)
        
    return client

def get_column(client, ssheet, wsheet, column, constraint={}):
    """Get the content of a specified column, optionally filtering on other columns"""
    
    # If the column specified is a name, find the corresponding index
    try:
        column = int(column)
    except ValueError:
        column = get_column_index(client,ssheet,wsheet,column)
    
    # Create a filter mask based on the supplied constraints
    filter = [True]*row_count(wsheet)
    for con_name, con_value in constraint.items():
        con_column = get_column(client,ssheet,wsheet,con_name)
        for i,value in enumerate(con_column):
            filter[i] &= (_to_unicode(value) == _to_unicode(con_value))
    
    # Get the content of the specified column index
    content_2d = get_cell_content(client, ssheet, wsheet, 0, column, 0, column)
    
    # Loop over the content and keep only the rows that have passed the constraint filters
    content = []
    for i,row in enumerate(content_2d):
        if filter[i]:
            content.append(row[0])
            
    return content
    
def get_column_index(client,ssheet,wsheet,name):
    """Get the index of the column with the specified name, or 0 if no column matches"""
    
    header = get_header(client,ssheet,wsheet)
    for i,column_name in enumerate(header):
        if _to_unicode(name) == _to_unicode(column_name):
            return (i+1)
    return 0

def get_header(client, ssheet, wsheet):
    """Return the column header of the supplied worksheet as a list"""
    
    header = get_row(client,ssheet,wsheet,1)
    return header

def get_key(object):
    """Get the unique gdocs key identifier for the supplied object"""
    return object.id.text.split('/')[-1]

def _get_query(title,exact_match):
    """Get a query object for the supplied parameters"""
    
    p = {}
    if title:
        # The urllib.quote method does not handle unicode, so encode the title in utf-8
        p['title'] = _from_unicode(title)
        if exact_match:
            p['title-exact'] = 'true'
        else:
            p['title-exact'] = 'false'
    
    return gdata.spreadsheet.service.DocumentQuery(params=p)

def get_row(client, ssheet, wsheet, row):
    """Get the content of a specified row index"""
    
    content = (get_cell_content(client, ssheet, wsheet, row, 0, row, 0) or [[]])
    return content[0]
    
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

def row_count(wsheet):
    """Get the number of rows in the worksheet"""
    return int(wsheet.row_count.text)

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
            client.UpdateCell(1,i+1,_to_unicode(header[i]),ss_key,ws_key)
    except:
        return False
    
    return True
        

    