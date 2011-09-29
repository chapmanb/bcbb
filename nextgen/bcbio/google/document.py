#!/usr/bin/python
"""Wrapper functions around Google Docs functionality for documents"""

import gdata.docs.service
import os
from gdata import MediaSource
from gdata import GDataEntry
import gdata.docs

from bcbio.google.connection import authenticate
from bcbio.google import (_from_unicode,_to_unicode)

def add_folder(client,folder_name,parent_folder=None):
    """Create a new folder"""
    # Check if the folder exists
    folder = get_folder(client,folder_name)
    if folder:
        return folder
    
    # Else, create it
    folder = client.CreateFolder(folder_name,parent_folder)
    return folder

def add_permission(client,doc,user,role_type="reader"):
    """Add the supplied permission to the document"""
    scope = gdata.docs.Scope(value=user, type='user')
    role = gdata.docs.Role(value=role_type)
    acl_entry = gdata.docs.DocumentListAclEntry(scope=scope, role=role)
    created_acl_entry = client.Post(acl_entry, doc.GetAclLink().href,converter=gdata.docs.DocumentListAclEntryFromString)
    
def add_spreadsheet(client,ssheet_title):
    """Create a new spreadsheet with the specified title"""
    new_entry = gdata.GDataEntry()
    new_entry.title = gdata.atom.Title(text=ssheet_title)
    category = client._MakeKindCategory(gdata.docs.service.SPREADSHEET_LABEL)
    new_entry.category.append(category)

    ssheet = client.Post(new_entry, '/feeds/documents/private/full')
    return ssheet

def get_client(encoded_credentials=None):
    """Get a DocsService client"""
    client = gdata.docs.service.DocsService()
    # If credentials were supplied, authenticate the client as well
    if encoded_credentials:
        authenticate(client,encoded_credentials)
        
    return client
  
def get_folder(client,folder_name):
    """Get a folder if it exists"""
    q = gdata.docs.service.DocumentQuery(categories=['folder'], params={'showfolders': 'true'})
    for entry in (client.Query(q.ToUri()).entry or []):
        if entry.title.text == folder_name:
            return entry
    return None
  
def move_to_folder(client,doc,folder):
    """Move a document into the supplied folder"""
    moved_doc = client.MoveIntoFolder(doc,folder)
    return moved_doc
    
def upload_file(client,fpath):
    """Upload a file to google docs
    """
    
    path, ext = os.path.splitext(fpath)
    ext = ext[1:]
    if not ext:
        return False
    ctype = gdata.docs.service.SUPPORTED_FILETYPES.get(ext.upper(),None)
    if not ctype:
        return False
    
    file_parent, file_name = os.path.split(path)
    ms = gdata.MediaSource(file_path=fpath,content_type=ctype)
    entry = client.Upload(ms, file_name)

    return entry
