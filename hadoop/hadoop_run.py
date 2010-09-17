#!/usr/bin/env python
"""Build up a hadoop job to process an input FASTA file.

Uses pydoop to manage interaction with the hadoop HDFS system.
"""
import os
import sys
import optparse
import subprocess
import contextlib

import pydoop.utils as pu
from pydoop.hdfs import hdfs

def main(script, in_dir, db_dir, out_dir):
    job_name = os.path.splitext(os.path.basename(script))[0]
    with _hdfs_filesystem() as (fs, lfs):
        script, in_dir, out_dir = setup_hdfs(job_name, script, in_dir,
                out_dir, fs, lfs)
        db_files_str, dbname = setup_db(job_name, db_dir, fs, lfs)
        hadoop_opts = {
          "mapred.job.name" : job_name,
          "hadoop.pipes.executable": script,
          "mapred.cache.files" : db_files_str,
          "mapred.create.symlink" : "yes",
          "hadoop.pipes.java.recordreader": "false",
          "hadoop.pipes.java.recordwriter": "true",
          "mapred.map.tasks": "2",
          "mapred.reduce.tasks": "2",
          "fasta.blastdb": dbname
          }

        cl = ["hadoop", "pipes"] + _cl_opts(hadoop_opts) + [
              "-program", script,
              "-input", in_dir, "-output", out_dir]
        subprocess.check_call(cl)
        process_output(fs, lfs, out_dir)

def process_output(fs, lfs, out_dir):
    for fname in (f['name'] for f in fs.list_directory(out_dir)):
        if os.path.basename(fname).startswith('part-'):
            in_handle = fs.open_file(fname)
            for line in in_handle:
                print line
            in_handle.close()

def setup_hdfs(work_dir, script, in_dir, out_dir, fs, lfs):
    """Add input output and script directories to hdfs.
    """
    _cleanup_workdir(fs, work_dir)
    out_info = []
    fs.create_directory(work_dir)
    # copy over input and files
    for lfile in [script, in_dir]:
        hdfs_ref = _hdfs_ref(work_dir, lfile)
        lfs.copy(lfile, fs, hdfs_ref)
        out_info.append(hdfs_ref)
    # create output directory
    hdfs_out = _hdfs_ref(work_dir, out_dir)
    out_info.append(hdfs_out)
    return out_info

def setup_db(work_dir, db_dir, fs, lfs):
    """Copy over BLAST database files, prepping them for map availability.
    """
    work_dir = os.path.join(work_dir, db_dir)
    fs.create_directory(work_dir)
    ref_info = []
    blast_dbs = []
    for fname in os.listdir(db_dir):
        blast_dbs.append(os.path.splitext(fname)[0])
        hdfs_ref = _hdfs_ref(work_dir, fname)
        lfs.copy(os.path.join(db_dir, fname), fs, hdfs_ref)
        ref_info.append("%s#%s" % (hdfs_ref, fname))
    return ",".join(ref_info), list(set(blast_dbs))[0]

@contextlib.contextmanager
def _hdfs_filesystem():
    """Retrieve references to the local and HDFS file system.

    Need to be able to specify host/port. For now, works off defaults.
    """
    fs = hdfs("default", 0)
    lfs = hdfs("", 0)
    try:
        yield fs, lfs
    finally:
        fs.close()
        lfs.close()

def _cleanup_workdir(fs, work_dir):
    if fs.exists(work_dir):
        fs.delete(work_dir)

def _hdfs_ref(work_dir, local_file):
    basename = os.path.basename(local_file)
    if basename == "":
        basename = os.path.basename(os.path.dirname(local_file))
    return os.path.join(work_dir, basename)

def _cl_opts(opts):
    cl = []
    for key, val in opts.iteritems():
        cl.append("-D")
        cl.append("%s=%s" % (key, val))
    return cl

if __name__ == "__main__":
    parser = optparse.OptionParser()
    opts, args = parser.parse_args()
    main(*args)

