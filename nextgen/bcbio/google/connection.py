#!/usr/bin/env python
"""Wrapper functions around the python gdata api for connecting and authenticating with the Google Docs service"""

import base64

def authenticate(client,credentials):
    
    login,pwd = _decode_credentials(credentials)
    if not login or not pwd:
        return False
    
    client.email = login
    client.password = pwd
    client.source = 'bcbb_nextgen_pipeline'
    client.ProgrammaticLogin()
    
    return True

def _decode_credentials(credentials):
    
    if not credentials:
        return None
    
    # Split the username and password
    return base64.b64decode(credentials).split(':',1);
 