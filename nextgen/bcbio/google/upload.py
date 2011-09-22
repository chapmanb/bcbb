#!/usr/bin/python
"""Upload a file to google docs
"""

import os
from gdata import MediaSource
import gdata.docs.service

def get_client():
    return gdata.docs.service.DocsService()

def upload_file(client,fpath):
    
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