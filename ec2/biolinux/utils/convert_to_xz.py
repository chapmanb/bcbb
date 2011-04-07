#!/usr/bin/env python
"""Convert gzipped files on s3 biodata to xz compression format.

This conversion is designed to save time and space for download.

Some download utilities to speed things up:
axel, aria2, lftp
"""
import os
import sys
import socket
import subprocess

import boto
import fabric.api as fabric

def main(bucket_name):
    conn = boto.connect_s3()
    bucket = conn.get_bucket("biodata")
    for s3_item in bucket.list("genomes/"):
        if s3_item.name.endswith(".gz"):
            print "xzipping", s3_item.name
            local_file = os.path.basename(s3_item.name)
            local_xz = "%s.xz" % os.path.splitext(local_file)[0]
            if not os.path.exists(local_xz):
                if not os.path.exists(local_file):
                    download_parallel(s3_item.generate_url(7200))
                    #s3_item.get_contents_to_filename(local_file)
                local_xz = gzip_to_xz(local_file)
            swap_s3_item(local_xz, bucket, s3_item)
            os.remove(local_xz)

def download_parallel(url):
    host = socket.gethostbyaddr(socket.gethostname())[0]
    user = os.environ["USER"]
    with fabric.settings(host_string="%s@%s" % (user, host)):
        ncores = fabric.run("cat /proc/cpuinfo | grep processor | wc -l")
        with fabric.cd(os.getcwd()):
            fabric.run("axel -a -n %s '%s'" % (ncores, url), shell=False)
            #fabric.run("aria2c -j %s -s %s '%s'" % (ncores, ncores, url),
            #           shell=False)

def swap_s3_item(xz_file, bucket, orig_s3_item):
    print " Uploading to S3"
    assert os.path.exists(xz_file)
    new_name = orig_s3_item.name.replace(".gz", ".xz")
    upload_script = os.path.join(os.path.dirname(__file__), "s3_multipart_upload.py")
    cl = ["python2.6", upload_script, xz_file, bucket.name, new_name]
    subprocess.check_call(cl)
    orig_s3_item.delete()

def gzip_to_xz(local_file):
    cl = ["gunzip", local_file]
    subprocess.check_call(cl)
    tar_file, _ = os.path.splitext(local_file)
    cl = ["xz", "-z", tar_file]
    subprocess.check_call(cl)
    return "%s.xz" % tar_file

if __name__ == "__main__":
    bucket_name = "biodata"
    main(bucket_name)
