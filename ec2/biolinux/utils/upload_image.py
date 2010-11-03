#!/usr/bin/env python
# Q&D EMI uploader

import os
import sys
import subprocess

prefix = "/scratch"
image = prefix+"/lucid-server-uec-amd64.img"

subprocess.Popen("euca-bundle-image -i "+image+ \
		" && euca-upload-bundle -b cloudbiolinux -m /tmp/"+os.path.basename(image)+".manifest.xml" \
		" && euca-register cloudbiolinux/"+os.path.basename(image)+".manifest.xml", shell=True)
