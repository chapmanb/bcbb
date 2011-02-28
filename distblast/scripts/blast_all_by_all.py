#!/usr/bin/env python
"""BLAST all proteins in a genome versus each other individually.

Prepares a matrix of identifies and BLAST scores for an all-versus-all
comparison of individual proteins versus each other.

Usage:
  blast_all_by_all.py <base_config.yaml> <org_config.yaml> <id_file.tsv>
"""
import os
import sys
import csv
import subprocess
import multiprocessing

import yaml
from scipy import sparse, io
from Bio import SeqIO

from bcbio.phylo import blast
from bcbio.picard.utils import chdir

def main(base_config_file, org_config_file, id_file):
    with open(base_config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    with open(id_file) as in_handle:
        ids = read_id_file(in_handle)
    data_out = "%s-all_search-data.tsv" % org_config["target_org"].replace(" ", "_")
    if not os.path.exists(data_out):
        with open(data_out, "w") as out_handle:
            prepare_data_file(out_handle, ids, org_config, config)
    write_out_matrices(ids, data_out)

def write_out_matrices(ids, data_out):
    base = os.path.splitext(data_out)[0].replace("-data", "")
    mat_file = "%s-scores.mat" % base
    with open(data_out) as in_handle:
        score_matrix, ident_matrix = get_matrices(in_handle, ids)
    io.savemat(mat_file, {"human_scores" : score_matrix,
                          "human_identities" : ident_matrix,
                          "human_ids" : ids})
    #id_file = "%s-ids.txt" % base
    #with open(id_file, "w") as out_handle:
    #    for i in ids:
    #        out_handle.write("%s\n" % i)

def get_matrices(in_handle, ids):
    pos_lookup = {}
    for pos, eid in enumerate(ids):
        pos_lookup[eid] = pos
    scores = sparse.lil_matrix((len(ids), len(ids)))
    idents = sparse.lil_matrix((len(ids), len(ids)))
    reader = csv.reader(in_handle, dialect="excel-tab")
    reader.next() # header
    for id1, id2, score, ident in reader:
        pos1 = pos_lookup[id1]
        pos2 = pos_lookup[id2]
        scores[pos1,pos2] = float(score)
        idents[pos1,pos2] = float(ident)
    return scores, idents

def prepare_data_file(out_handle, ids, org_config, config):
    writer = csv.writer(out_handle, dialect="excel-tab")
    seq_recs = SeqIO.index(org_config["search_file"], "fasta")
    search_db = make_search_db(seq_recs, ids, org_config["target_org"], config["work_dir"])
    writer.writerow(["rec1", "rec2", "score", "identity"])
    pool = multiprocessing.Pool(int(config["num_cores"]))
    results = pool.imap(blast_seqs,
                        ((i, seq_recs[i], search_db, config) for i in ids))
    for info in results:
        for id1, id2, identity, score in info:
            writer.writerow([id1, id2, score, identity])

def make_search_db(seq_recs, ids, target_org, tmp_dir):
    search_db = "%s-db.fa" % target_org.replace(" ", "_")
    db_name = os.path.splitext(search_db)[0]
    with chdir(tmp_dir):
        with open(search_db, "w") as out_handle:
            SeqIO.write((seq_recs[i] for i in ids), out_handle, "fasta")
        cl = ["makeblastdb", "-in", search_db,
              "-dbtype", "prot",
              "-out", db_name,
              "-title", target_org]
        subprocess.check_call(cl)
    return os.path.join(tmp_dir, db_name)

def blast_seqs(args):
    """Blast two sequences, returning the score and identity.
    """
    def do_work(rec_id, rec, org_db, config):
        return blast.blast_hit_list(rec_id, rec, org_db, config["evalue_thresh"], config["work_dir"])
    return do_work(*args)

def read_id_file(in_handle):
    in_handle.next() # header
    return [l.split("\t")[0] for l in in_handle]

if __name__ == "__main__":
    main(*sys.argv[1:])
