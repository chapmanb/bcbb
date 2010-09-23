#!/usr/bin/env python
"""Build up a hadoop job to process an input FASTA file.

Uses pydoop to manage interaction with the hadoop HDFS system.
"""
import os
import sys
import csv
import glob
import json
import shutil
import optparse
import subprocess
import contextlib

import yaml
import pydoop.utils as pu
from pydoop.hdfs import hdfs

from bcbio.phylo import blast

def main(script, org_config_file, config_file, in_dir, local_out_dir):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    if os.path.exists(in_dir):
        shutil.rmtree(in_dir)
    os.makedirs(in_dir)
    shutil.copy(org_config["search_file"], in_dir)

    job_name = os.path.splitext(os.path.basename(script))[0]
    with _hdfs_filesystem() as (fs, lfs):
        script, in_dir, out_dir = setup_hdfs(job_name, script, in_dir,
                local_out_dir, fs, lfs)
        db_files_str, dbnames, org_names = setup_db(job_name,
                config['db_dir'], org_config['target_org'], fs, lfs)
        hadoop_opts = {
          "mapred.job.name" : job_name,
          "hadoop.pipes.executable": script,
          "mapred.cache.files" : db_files_str,
          "mapred.create.symlink" : "yes",
          "hadoop.pipes.java.recordreader": "false",
          "hadoop.pipes.java.recordwriter": "true",
          "mapred.map.tasks": "2",
          "mapred.reduce.tasks": "2",
          "fasta.blastdb": ",".join(dbnames)
          }

        cl = ["hadoop", "pipes"] + _cl_opts(hadoop_opts) + [
              "-program", script,
              "-input", in_dir, "-output", out_dir]
        subprocess.check_call(cl)
        process_output(fs, lfs, out_dir, local_out_dir,
                org_config['target_org'], org_names)

def process_output(fs, lfs, out_dir, local_out_dir, target_org, org_names):
    """Convert JSON output into tab delimited score and ID files.
    """
    # setup the output files
    if not os.path.exists(local_out_dir):
        os.makedirs(local_out_dir)
    base = target_org.replace(" ", "_")
    id_file = os.path.join(local_out_dir, "%s-ids.tsv" % base)
    score_file = os.path.join(local_out_dir, "%s-scores.tsv" % base)
    with open(id_file, "w") as id_handle:
        with open(score_file, "w") as score_handle:
            id_writer = csv.writer(id_handle, dialect="excel-tab")
            score_writer = csv.writer(score_handle, dialect="excel-tab")
            header = [""] + org_names
            id_writer.writerow(header)
            score_writer.writerow(header)
            # loop through the hadoop files and reformat each
            # into tab delimited output
            for fname in (f['name'] for f in fs.list_directory(out_dir)):
                if os.path.basename(fname).startswith('part-'):
                    in_handle = fs.open_file(fname)
                    for line in in_handle:
                        cur_id, info = line.split("\t")
                        data = json.loads(info)
                        id_writer.writerow([cur_id] + data['ids'])
                        score_writer.writerow([cur_id] + data['scores'])
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

def setup_db(work_dir, db_dir, target_org, fs, lfs):
    """Copy over BLAST database files, prepping them for map availability.
    """
    (_, org_names, db_refs) = blast.get_org_dbs(db_dir, target_org)
    work_dir = os.path.join(work_dir, db_dir)
    fs.create_directory(work_dir)
    db_refs = db_refs
    ref_info = []
    blast_dbs = []
    for db_path in db_refs:
        blast_dbs.append(os.path.basename(db_path))
        for fname in glob.glob(db_path + ".[p|n]*"):
            hdfs_ref = _hdfs_ref(work_dir, fname)
            lfs.copy(fname, fs, hdfs_ref)
            ref_info.append("%s#%s" % (hdfs_ref, os.path.basename(hdfs_ref)))
    return ",".join(ref_info), blast_dbs, org_names

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

