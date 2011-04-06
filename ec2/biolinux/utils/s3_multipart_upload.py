#!/usr/bin/env python
"""Split large file into multiple pieces for upload to S3.

S3 only supports 5Gb files for uploading directly, so for larger CloudBioLinux
box images we need to use boto's multipart file support.

This parallelizes the task over available cores using multiprocessing.

Usage:
  s3_multipart_upload.py <file_to_transfer> <bucket_name> [<s3_key_name>]
    if <s3_key_name> is not specified, the filename will be used.

    --norr -- Do not use reduced redundancy storage.
    --nopublic -- Do not make files public.

    Files are stored at cheaper reduced redundancy storage and made public by
    default.
"""
import os
import sys
import glob
import subprocess
import contextlib
import functools
import multiprocessing
from multiprocessing.pool import IMapIterator
from optparse import OptionParser

import boto

def main(transfer_file, bucket_name, s3_key_name=None, use_rr=True,
         make_public=True):
    if s3_key_name is None:
        s3_key_name = os.path.basename(transfer_file)
    conn = boto.connect_s3()
    bucket = conn.lookup(bucket_name)
    mb_size = os.path.getsize(transfer_file) / 1e6
    if mb_size < 60:
        _standard_transfer(bucket, s3_key_name, transfer_file, use_rr)
    else:
        _multipart_upload(bucket, s3_key_name, transfer_file, mb_size, use_rr)
    s3_key = bucket.get_key(s3_key_name)
    if make_public:
        s3_key.set_acl("public-read")

def upload_cb(complete, total):
    sys.stdout.write(".")
    sys.stdout.flush()

def _standard_transfer(bucket, s3_key_name, transfer_file, use_rr):
    print " Upload with standard transfer, not multipart",
    new_s3_item = bucket.new_key(s3_key_name)
    new_s3_item.set_contents_from_filename(transfer_file, reduced_redundancy=use_rr,
                                           cb=upload, num_cb=10)
    print

def map_wrap(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return apply(f, *args, **kwargs)
    return wrapper

def mp_from_ids(mp_id, mp_keyname, mp_bucketname):
    """Get the multipart upload from the bucket and multipart IDs.

    This allows us to reconstitute a connection to the upload
    from within multiprocessing functions.
    """
    conn = boto.connect_s3()
    bucket = conn.lookup(mp_bucketname)
    mp = boto.s3.multipart.MultiPartUpload(bucket)
    mp.key_name = mp_keyname
    mp.id = mp_id
    return mp

@map_wrap
def transfer_part(mp_id, mp_keyname, mp_bucketname, i, part):
    """Transfer a part of a multipart upload. Designed to be run in parallel.
    """
    mp = mp_from_ids(mp_id, mp_keyname, mp_bucketname)
    print " Transferring", i, part
    with open(part) as t_handle:
        mp.upload_part_from_file(t_handle, i+1)
    os.remove(part)

def _multipart_upload(bucket, s3_key_name, tarball, mb_size, use_rr=True):
    """Upload large files using Amazon's multipart upload functionality.
    """
    def split_file(in_file, mb_size):
        prefix = os.path.join(os.path.dirname(in_file),
                              "%sS3PART" % (os.path.basename(s3_key_name)))
        split_size = int(min(mb_size / 10.0, 250))
        if not os.path.exists("%saa" % prefix):
            cl = ["split", "-b%sm" % split_size, in_file, prefix]
            subprocess.check_call(cl)
        return sorted(glob.glob("%s*" % prefix))

    mp = bucket.initiate_multipart_upload(s3_key_name, reduced_redundancy=use_rr)
    with multimap() as pmap:
        for _ in pmap(transfer_part, ((mp.id, mp.key_name, mp.bucket_name, i, part) for (i, part) in
                                      enumerate(split_file(tarball, mb_size)))):
            pass
    mp.complete_upload()

@contextlib.contextmanager
def multimap(cores=None):
    """Provide multiprocessing imap like function.

    The context manager handles setting up the pool, worked around interrupt issues
    and terminating the pool on completion.
    """
    if cores is None:
        cores = max(multiprocessing.cpu_count() - 1, 1)
    def wrapper(func):
        def wrap(self, timeout=None):
            return func(self, timeout=timeout if timeout is not None else 1e100)
        return wrap
    IMapIterator.next = wrapper(IMapIterator.next)
    pool = multiprocessing.Pool(cores)
    yield pool.imap
    pool.terminate()

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-r", "--norr", dest="use_rr",
                      action="store_false", default=True)
    parser.add_option("-p", "--nopublic", dest="make_public",
                      action="store_false", default=True)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    kwags = dict(use_rr=options.use_rr, make_public=options.make_public)
    main(*args, **kwargs)
