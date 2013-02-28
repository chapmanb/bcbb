#!/usr/bin/env python
"""Perform genome wide BLAST comparisons of an organism against other genomes.

Usage:
    blast_cross_orgs.py <organism config> <YAML config>

This requires a set of BLAST databases setup by 'retrieve_org_dbs.py'.

Requires:
    - NCBI's blast+
    - Biopython libraries
"""
import os
import sys
import csv
import itertools
import multiprocessing

import yaml
from Bio import SeqIO

from bcbio.phylo import blast

def main(org_config_file, config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    if not os.path.exists(config['work_dir']):
        os.makedirs(config['work_dir'])
    (_, org_names, db_refs) = blast.get_org_dbs(config['db_dir'],
            org_config['target_org'])
    id_file, score_file = setup_output_files(org_config['target_org'])
    with open(org_config['search_file']) as in_handle:
        with open(id_file, "w") as id_out_handle:
            with open(score_file, "w") as score_out_handle:
                id_writer = csv.writer(id_out_handle, dialect='excel-tab')
                score_writer = csv.writer(score_out_handle, dialect='excel-tab')
                header = [""] + org_names
                id_writer.writerow(header)
                score_writer.writerow(header)
                _do_work(db_refs, in_handle, id_writer, score_writer, config)

def _do_work(db_refs, in_handle, id_writer, score_writer, config):
    cores = int(config["num_cores"])
    pool = multiprocessing.Pool(cores)
    for rec_group in partition_all(cores * 100, SeqIO.parse(in_handle, "fasta")):
        for out in pool.imap(_process_wrapper,
                             ((rec, db_refs, config['work_dir'], config.get("blast_cmd"))
                              for rec in rec_group)):
            id_writer.writerow([out["cur_id"]] + out["cmp_id"])
            score_writer.writerow([out["cur_id"]] + out["cmp_score"])

def partition_all(n, iterable):
    """Split into lazy chunks of n size.
    http://stackoverflow.com/questions/5129102/python-equivalent-to-clojures-partition-all
    """
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, n))
        if not chunk:
            break
        yield chunk

def _process_wrapper(args):
    try:
        return process_blast(*args)
    except KeyboardInterrupt:
        raise Exception

def process_blast(rec, db_refs, tmp_dir, blast_cmd):
    """Run a BLAST writing results to shared files.
    """
    cur_id, id_info, score_info = blast.blast_top_hits(rec.id, rec.format("fasta"),
            db_refs, tmp_dir, blast_cmd)
    print cur_id
    return {"cmp_id": id_info,
            "cmp_score": score_info,
            "cur_id": cur_id}

def setup_output_files(target_org):
    base = target_org.replace(" ", "_")
    id_file = "%s-ids.tsv" % base
    eval_file = "%s-scores.tsv" % base
    return id_file, eval_file

if __name__ == "__main__":
    main(*sys.argv[1:])
