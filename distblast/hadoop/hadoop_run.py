#!/usr/bin/env python
"""Build up a hadoop job to process an input FASTA file.

This handles copying over all of the input files, scripts and databases.
After running the Hadoop job, the output files are retrieved and reformatted
to be identical to a local run.
"""
import os
import sys
import csv
import glob
import json
import shutil
import optparse
import subprocess

import yaml

from bcbio.phylo import blast

def main(script, org_config_file, config_file, in_dir, out_dir):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    if os.path.exists(in_dir):
        shutil.rmtree(in_dir)
    os.makedirs(in_dir)
    shutil.copy(org_config["search_file"], in_dir)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    job_name = os.path.splitext(os.path.basename(script))[0]
    print "Creating working directories in HDFS"
    hdfs_script, hdfs_in_dir, hdfs_out_dir = setup_hdfs(job_name, script,
                                                        in_dir, out_dir)
    print "Copying organism database files to HDFS"
    db_files, dbnames, org_names = setup_db(job_name, config['db_dir'],
                                            org_config['target_org'])
    print "Running script on Hadoop"
    if script.endswith("streaming.py"):
        run_hadoop_streaming(job_name, script, hdfs_in_dir, hdfs_out_dir,
                             db_files, dbnames)
    else:
        run_hadoop(job_name, hdfs_script, hdfs_in_dir, hdfs_out_dir,
                   db_files, dbnames)
    print "Processing output files"
    process_output(hdfs_out_dir, out_dir, org_config['target_org'],
                   org_names)

def run_hadoop_streaming(job_name, script, hdfs_in_dir, hdfs_out_dir, db_files, dbnames):
    """Run a hadoop streaming job with dumbo.
    """
    cachefiles = []
    for db_file in db_files:
        cachefiles.extend(("-cachefile", db_file))
    hadoop = os.path.join(os.environ["HADOOP_HOME"], "bin", "hadoop")
    hadoop = os.environ["HADOOP_HOME"]
    cl = ["dumbo", script, "-hadoop", hadoop, "-input", hdfs_in_dir,
          "-output", hdfs_out_dir, "-name", job_name,
          "-inputformat", "com.lifetech.hadoop.streaming.FastaInputFormat",
          "-libjar", os.path.abspath("jar/bioseq-0.0.1.jar"),
          "-outputformat", "text",
          "-cmdenv", "fasta_blastdb=%s" % ','.join(dbnames)] + cachefiles
    subprocess.check_call(cl)

def run_hadoop_streaming_mrjob(job_name, script, in_dir, db_files_str, dbnames):
    jobconf_opts = {
        "mapred.job.name" : job_name,
        "mapred.cache.files" : db_files_str,
        "mapred.create.symlink" : "yes",
        "fasta.blastdb": ",".join(dbnames)
        #"mapred.task.timeout": "60000", # useful for debugging
        }
    hadoop_opts = {
        "-inputformat": "com.lifetech.hadoop.streaming.FastaInputFormat",
        "-libjars" : os.path.abspath("jars/bioseq-0.0.1.jar")
        }
    def _mrjob_opts(opts, type, connect):
        out = []
        for k, v in opts.iteritems():
            if connect == " ":
                val = "'%s %s'" % (k, v)
            else:
                val = "%s%s%s" % (k, connect, v)
            out.append("%s=%s" % (type, val))
        return out
    cl = [sys.executable, script, "-r", "hadoop"] + \
         _mrjob_opts(jobconf_opts, "--jobconf", "=") + \
         _mrjob_opts(hadoop_opts, "--hadoop-arg", " ") + \
         glob.glob(os.path.join(in_dir, "*"))
    subprocess.check_call(cl)

def run_hadoop(job_name, hdfs_script, hdfs_in_dir, hdfs_out_dir,
               db_files, dbnames):
    hadoop_opts = {
      "mapred.job.name" : job_name,
      "hadoop.pipes.executable": hdfs_script,
      "mapred.cache.files" : ",".join(db_files),
      "mapred.create.symlink" : "yes",
      "hadoop.pipes.java.recordreader": "false",
      "hadoop.pipes.java.recordwriter": "true",
      "mapred.map.tasks": "2",
      "mapred.reduce.tasks": "2",
      #"mapred.task.timeout": "60000", # useful for debugging
      "fasta.blastdb": ",".join(dbnames)
      }
    cl = ["hadoop", "pipes"] + _cl_opts(hadoop_opts) + [
          "-program", hdfs_script,
          "-input", hdfs_in_dir, "-output", hdfs_out_dir]
    subprocess.check_call(cl)

def _read_hadoop_out(out_dir, work_dir):
    """Reformat Hadoop output files for tab delimited output.
    """
    cl = ["hadoop", "fs", "-ls", os.path.join(out_dir, "part-*")]
    p = subprocess.Popen(cl, stdout=subprocess.PIPE)
    (out, _) = p.communicate()
    for fname in sorted([l.split()[-1] for l in out.split("\n")
                         if l and not l.startswith("Found")]):
        local_fname = os.path.join(work_dir, os.path.basename(fname))
        cl = ["hadoop", "fs", "-get", fname, local_fname]
        subprocess.check_call(cl)
        with open(local_fname) as in_handle:
            for line in in_handle:
                cur_id, info = line.split("\t")
                data = json.loads(info)
                data["cur_id"] = cur_id
                yield data
        os.remove(local_fname)

def process_output(hdfs_out_dir, local_out_dir, target_org, org_names):
    """Convert JSON output into tab delimited score and ID files.
    """
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
            for data in _read_hadoop_out(hdfs_out_dir, local_out_dir):
                id_writer.writerow([data["cur_id"]] + data['ids'])
                score_writer.writerow([data["cur_id"]] + data['scores'])

def setup_hdfs(hdfs_work_dir, script, in_dir, out_dir):
    """Add input, output and script directories to hdfs.
    """
    cl = ["hadoop", "fs", "-test", "-d", hdfs_work_dir]
    result = subprocess.call(cl)
    if result == 0:
        cl = ["hadoop", "fs", "-rmr", hdfs_work_dir]
        subprocess.check_call(cl)
    cl = ["hadoop", "fs", "-mkdir", hdfs_work_dir]
    subprocess.check_call(cl)
    out_info = []
    # copy over input and files
    for lfile in [script, in_dir]:
        hdfs_ref = _hdfs_ref(hdfs_work_dir, lfile)
        cl = ["hadoop", "fs", "-put", lfile, hdfs_ref]
        subprocess.check_call(cl)
        out_info.append(hdfs_ref)
    hdfs_out = _hdfs_ref(hdfs_work_dir, out_dir)
    out_info.append(hdfs_out)
    return out_info

def setup_db(work_dir_base, db_dir, target_org):
    """Copy over BLAST database files, prepping them for map availability.
    """
    (_, org_names, db_refs) = blast.get_org_dbs(db_dir, target_org)
    work_dir = os.path.join(work_dir_base, db_dir)
    cl = ["hadoop", "fs", "-mkdir", work_dir]
    subprocess.check_call(cl)
    ref_info = []
    blast_dbs = []
    for db_path in db_refs:
        blast_dbs.append(os.path.basename(db_path))
        for fname in glob.glob(db_path + ".[p|n]*"):
            hdfs_ref = _hdfs_ref(work_dir, fname)
            cl = ["hadoop", "fs", "-put", fname, hdfs_ref]
            subprocess.check_call(cl)
            ref_info.append("%s#%s" % (hdfs_ref, os.path.basename(hdfs_ref)))
    return ref_info, blast_dbs, org_names

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
