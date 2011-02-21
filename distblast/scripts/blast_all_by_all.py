#!/usr/bin/env python
"""BLAST all proteins in a genome versus each other individually.

Prepares a matrix of identifies and BLAST scores for an all-versus-all
comparison of individual proteins versus each other.

Usage:
  blast_all_by_all.py <base_config.yaml> <org_config.yaml> <id_file.tsv>
"""
import sys
import csv
import subprocess
import multiprocessing

import yaml
from Bio import SeqIO

from bcbio.phylo import blast

def main(base_config_file, org_config_file, id_file):
    with open(base_config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    with open(id_file) as in_handle:
        ids = read_id_file(in_handle)
    seq_recs = SeqIO.index(org_config["search_file"], "fasta")
    out = OutputManager(org_config["target_org"])
    out.write_header(ids)
    for i, id1 in enumerate(ids):
        compares = [""] * (i) + ids[i:]
        pool = multiprocessing.Pool(int(config["num_cores"]))
        assert len(compares) == len(ids)
        results = pool.imap(blast_seqs,
                            ((seq_recs[id1], seq_recs.get(id2, ""), config) for id2 in compares))
        out.write_id(id1, list(results))

class OutputManager:
    def __init__(self, org_name):
        safe_org_name = org_name.replace(" ", "_")
        score_out = open("%s-all_by_all-scores.tsv" % safe_org_name, "w")
        self.score_writer = csv.writer(score_out, dialect="excel-tab")
        ident_out = open("%s-all_by_all-identities.tsv" % safe_org_name, "w")
        self.ident_writer = csv.writer(ident_out, dialect="excel-tab")

    def write_header(self, ids):
        self.score_writer.writerow([""] + ids)
        self.ident_writer.writerow([""] + ids)

    def write_id(self, name, results):
        self.score_writer.writerow([name] + [x[0] for x in results])
        self.ident_writer.writerow([name] + [x[1] for x in results])

def blast_seqs(args):
    """Blast two sequences, returning the score and identity.
    """
    def do_work(rec1, rec2, config):
        if rec2:
            identity, score = blast.blast_two_seqs(rec1, rec2, config["work_dir"])
        else:
            score, identity = "", ""
        return score, identity
    return do_work(*args)

def read_id_file(in_handle):
    in_handle.next() # header
    return [l.split("\t")[0] for l in in_handle]

if __name__ == "__main__":
    main(*sys.argv[1:])
