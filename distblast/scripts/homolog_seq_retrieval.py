#!/usr/bin/env python
"""Retrieve sequences for a group of homologs based on a gene list.

Usage:
    homolog_seq_retrieval.py <base config> <org config> <id result file>
"""
import sys
import csv
import os
import re

import yaml
from Bio import SeqIO

def main(base_file, org_file, result_file):
    with open(base_file) as in_handle:
        base_config = yaml.load(in_handle)
    with open(org_file) as in_handle:
        org_config = yaml.load(in_handle)
    out_dir = os.path.join(os.getcwd(), "%s-seqs" %
            org_config['target_org'].replace(" ", "-"))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    id_list = read_gene_list(org_config['gene_list'])
    id_list = _build_alt_tx_res(id_list)
    with open(result_file) as in_handle:
        reader = csv.reader(in_handle, dialect='excel-tab')
        org_list = reader.next()[1:]
        id_indexes = prepare_indexes(org_config['search_file'],
                base_config['db_dir'], org_list)
        print "Retrieving sequences from list"
        all_matches = dict()
        for cur_id, org_ids in ((l[0], l[1:]) for l in reader):
            match_index = _transcript_matches(cur_id, id_list)
            if match_index > -1:
                print cur_id
                out_file = os.path.join(out_dir, "%s.fa" %
                        cur_id.replace(".", "_"))
                output_seqs([cur_id] + org_ids,
                    [org_config['target_org']] + org_list,
                    id_indexes, out_file)
                all_matches[match_index] = ""
    print "Not found", [id_list[i][0] for i in set(range(len(id_list))) -
            set(all_matches.keys())]

def output_seqs(ids, orgs, indexes, out_file):
    """Output the sequences we are interested in retrieving.
    """
    with open(out_file, "w") as out_handle:
        SeqIO.write(_all_seqs(ids, orgs, indexes), out_handle, "fasta")

def _all_seqs(ids, orgs, indexes):
    """Lazy generator of sequences from our indexes, with IDs properly set.
    """
    for i, cur_id in enumerate(ids):
        if cur_id:
            rec = indexes[i][cur_id]
            rec.id = cur_id
            rec.description = orgs[i]
            yield rec

def _build_alt_tx_res(id_list):
    """Create regular expressions to find alternative transcripts.

    Matches alternative transcripts of base.1 of the form:
        base.1a
        base.1.1
        base.1a.1
    """
    final_list = []
    for cur_id in id_list:
        match_re = re.compile(r"^" + cur_id + r"([a-z]+|\.\d+)(\.\d+)?$")
        final_list.append((cur_id, match_re))
    return final_list

def _transcript_matches(cur_id, id_list):
    """Check if a transcript matches one of the alternatively spliced names.
    """
    for index, (test_id, test_re) in enumerate(id_list):
        if cur_id == test_id or test_re.search(cur_id):
            return index
    return  -1

def prepare_indexes(base_file, db_dir, orgs):
    """Prepare easy to retrieve Biopython indexes from Fasta inputs.
    """
    print "Preparing fasta indexes"
    indexes = []
    for fname in [base_file] + [_get_org_fasta(db_dir, o) for o in orgs]:
        normalizer = IdNormalizer()
        indexes.append(SeqIO.index(fname, "fasta",
            key_function=normalizer.normalize))
        normalizer.finished = True
    return indexes

def _get_org_fasta(db_dir, org):
    """Retrieve the FASTA file for an organism database.
    """
    base_file = os.path.join(db_dir, "organism_dbs.txt")
    with open(base_file) as in_handle:
        for line in in_handle:
            if line.startswith(org):
                parts = line.strip().split()
                return os.path.join(db_dir, "%s.fa" % parts[-1])
    raise ValueError("Did not find fasta file for %s" % org)

def read_gene_list(in_file):
    """Read list of genes of interest.
    """
    genes = []
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for parts in reader:
            tid = parts[-1]
            genes.append(tid.strip())
    return genes

class IdNormalizer:
    def __init__(self):
        self._seen_ids = {}
        self._index = 0
        self.finished = False

    def normalize(self, id_info):
        if id_info.startswith("gi|"):
            parts = [p for p in id_info.split("|") if p]
            id_info = parts[-1]
        if not self.finished:
            try:
                self._seen_ids[id_info]
                self._index += 1
                return self._index
            except KeyError:
                self._seen_ids[id_info] = ""
                return id_info
        return id_info

if __name__ == "__main__":
    apply(main, sys.argv[1:])
