#!/usr/bin/env python
"""Split large file into multiple pieces for upload to S3.

S3 only supports 5Gb files for uploading directly, so for larger CloudBioLinux
box images we need to use boto's multipart file support.
"""
import os
import sys
import glob
import subprocess

import boto

def main(transfer_file, bucket_name):
    conn = boto.connect_s3()
    bucket = conn.lookup(bucket_name)
    mp = bucket.initiate_multipart_upload(transfer_file)

    for i, part in enumerate(split_file(transfer_file)):
        print "Transferring", part
        with open(part) as t_handle:
            mp.upload_part_from_file(t_handle, i+1)
        os.remove(part)
    mp.complete_upload()

def split_file(in_file):
    prefix = "S3PART"
    cl = ["split", "-b100m", in_file, prefix]
    subprocess.check_call(cl)
    return sorted(glob.glob("%s*" % prefix))

if __name__ == "__main__":
    main(*sys.argv[1:])
