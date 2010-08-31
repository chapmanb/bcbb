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
import subprocess
import contextlib
import multiprocessing
import xml.parsers.expat

import yaml
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

fupdate_lock = multiprocessing.Lock()

def main(org_config_file, config_file):
    # fasta_ref
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(org_config_file) as in_handle:
        org_config = yaml.load(in_handle)
    if not os.path.exists(config['work_dir']):
        os.makedirs(config['work_dir'])
    (_, db_refs) = get_org_dbs(config['db_dir'], org_config['target_org'])
    id_file, score_file = setup_output_files(org_config['target_org'],
            [r[0] for r in db_refs])
    file_info = [id_file, score_file]
    pool = multiprocessing.Pool(int(config['num_cores']))
    with open(org_config['search_file']) as in_handle:
        pool.map(_process_wrapper, ((config, org_config, i, rec, db_refs, file_info)
            for (i, rec) in enumerate(SeqIO.parse(in_handle, "fasta"))))

def process_blast(config, org_config, i, rec, db_refs, file_info):
    cur_id = _normalize_id(rec.id)
    id_info = [cur_id]
    score_info = [cur_id]
    for xorg, xdb in db_refs:
        with _blast_filenames(config, org_config, xorg, i) as (in_file, out_file):
            SeqIO.write([rec], in_file, "fasta")
            out_id, out_eval = compare_by_blast(in_file, xdb, out_file)
            id_info.append(out_id)
            score_info.append(out_eval)
    with fupdate_lock:
        id_file, score_file = file_info
        for fname, fvals in [(id_file, id_info), (score_file, score_info)]:
            with open(fname, "a") as out_handle:
                writer = csv.writer(out_handle, dialect='excel-tab')
                writer.writerow(fvals)
        print cur_id

def _process_wrapper(args):
    try:
        return process_blast(*args)
    except KeyboardInterrupt:
        raise Exception

def compare_by_blast(input_ref, xref_db, blast_out):
    """Compare all genes in an input file to the output database.
    """
    cl = NcbiblastpCommandline(query=input_ref, db=xref_db, out=blast_out,
            outfmt=5, num_descriptions=1, num_alignments=0)
    subprocess.check_call(str(cl).split())
    with open(blast_out) as blast_handle:
        try:
            rec = NCBIXML.read(blast_handle)
        except xml.parsers.expat.ExpatError:
            rec = None
        if rec and len(rec.descriptions) > 0:
            id_info = _normalize_id(rec.descriptions[0].title.split()[1])
            return id_info, rec.descriptions[0].bits
        else:
            return "", 0

def _normalize_id(id_info):
    if id_info.startswith("gi|"):
        parts = [p for p in id_info.split("|") if p]
        id_info = parts[-1]
    return id_info

def get_org_dbs(db_dir, target_org):
    """Retrieve references to fasta and BLAST databases for included organisms.
    """
    fasta_ref = None
    db_refs = []
    with open(os.path.join(db_dir, "organism_dbs.txt")) as org_handle:
        for line in org_handle:
            org, db = line.rstrip("\r\n").split("\t")
            if db:
                if org == target_org:
                    assert fasta_ref is None
                    fasta_ref = os.path.join(db_dir, "%s.fa" % db)
                db_refs.append((org, os.path.join(db_dir, db)))
    assert fasta_ref is not None, "Did not find base organism database"
    return fasta_ref, db_refs

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

@contextlib.contextmanager
def _blast_filenames(config, org_config, xorg, i):
    """Create files needed for blast and cleanup on closing.
    """
    base_name = "%s_%s_%s" % (org_config['target_org'].replace(" ", "_"),
                              xorg.replace(" ", "_"), i)
    in_file = os.path.join(config['work_dir'], "%s-in.fa" % base_name)
    out_file = os.path.join(config['work_dir'], "%s-out.blast" % base_name)
    try:
        yield (in_file, out_file)
    finally:
        for fname in [in_file, out_file]:
            if os.path.exists(fname):
                os.remove(fname)

@contextlib.contextmanager
def _chdir(new_dir):
    orig_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(orig_dir)

if __name__ == "__main__":
    main(*sys.argv[1:])
