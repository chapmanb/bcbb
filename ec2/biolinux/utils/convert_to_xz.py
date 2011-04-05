#!/usr/bin/env python
"""Convert gzipped files on s3 biodata to xz compression format.

This conversion is designed to save time and space for download.
"""
import os
import sys
import subprocess

import boto

def main(bucket_name):
    conn = boto.connect_s3()
    bucket = conn.get_bucket("biodata")
    for s3_item in bucket.list("genomes/"):
        if s3_item.name.endswith(".gz"):
            print "xzipping", s3_item.name
            local_file = os.path.basename(s3_item.name)
            if not os.path.exists(local_file):
                s3_item.get_contents_to_filename(local_file)
            local_xz = gzip_to_xz(local_file)
            swap_s3_item(local_xz, bucket, s3_item)
            os.remove(local_xz)

def upload_cb(complete, total):
    sys.stdout.write(".")
    sys.stdout.flush()

def swap_s3_item(xz_file, bucket, orig_s3_item):
    print " Uploading to S3",
    assert os.path.exists(xz_file)
    new_name = orig_s3_item.name.replace(".gz", ".xz")
    new_s3_item = bucket.new_key(new_name)
    new_s3_item.set_contents_from_filename(xz_file, reduced_redundancy=True,
                                           cb=upload_cb, num_cb=10)
    new_s3_item.make_public()
    orig_s3_item.delete()
    print

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
