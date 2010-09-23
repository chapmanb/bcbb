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
import multiprocessing

import yaml
from Bio import SeqIO

from bcbio.phylo import blast

fupdate_lock = multiprocessing.Lock()

def main(org_config_file, config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    if not os.path.exists(config['work_dir']):
        os.makedirs(config['work_dir'])
    (_, org_names, db_refs) = blast.get_org_dbs(config['db_dir'],
            org_config['target_org'])
    id_file, score_file = setup_output_files(org_config['target_org'],
            org_names)
    file_info = [id_file, score_file]
    pool = multiprocessing.Pool(int(config['num_cores']))
    with open(org_config['search_file']) as in_handle:
        pool.map(_process_wrapper,
                ((rec, db_refs, file_info, config['work_dir'])
                    for rec in SeqIO.parse(in_handle, "fasta")))

def _process_wrapper(args):
    try:
        return process_blast(*args)
    except KeyboardInterrupt:
        raise Exception

def process_blast(rec, db_refs, file_info, tmp_dir):
    """Run a BLAST writing results to shared files.
    """
    cur_id, id_info, score_info = blast.blast_top_hits(rec.id, rec.format("fasta"),
            db_refs, tmp_dir)
    with fupdate_lock:
        id_file, score_file = file_info
        for fname, fvals in [(id_file, id_info), (score_file, score_info)]:
            with open(fname, "a") as out_handle:
                writer = csv.writer(out_handle, dialect='excel-tab')
                writer.writerow([cur_id] + fvals)
        print cur_id

def setup_output_files(target_org, cmp_orgs):
    base = target_org.replace(" ", "_")
    id_file = "%s-ids.tsv" % base
    eval_file = "%s-scores.tsv" % base
    for fname in [id_file, eval_file]:
        with open(fname, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            header = [""] + cmp_orgs
            writer.writerow(header)
    return id_file, eval_file

if __name__ == "__main__":
    main(*sys.argv[1:])
