#!/usr/bin/env python
"""Filter an output file, removing alternative transcripts based on names.

Filter ID and score output files from a distributed BLAST, returning
the longest transcript for a gene as an evolutionary transcript sample for
that gene.

Usage:
    filter_by_transcript.py <org config file>
"""
import sys
import os
import csv
import re
import collections
from optparse import OptionParser

import yaml
from Bio import SeqIO

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    target_files = files_to_filter(config["target_org"])
    names_to_include = get_representative_txs(config['search_file'])
    for f in target_files:
        filter_file(f, names_to_include)

def filter_file(to_filter, names_to_include):
    new_ext = "longtxs"
    base, ext = os.path.splitext(to_filter)
    out_file = "%s-%s%s" % (base, new_ext, ext)
    with open(to_filter) as in_handle:
        with open(out_file, "w") as out_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(reader.next())
            for parts in reader:
                if parts[0] in names_to_include:
                    writer.writerow(parts)

def get_representative_txs(in_file):
    """Retrieve the largest transcript for each gene, using Ensembl headers.

    This relies on Ensembl header structure in FASTA files. For all genes,
    the largest of potentially many alternative transcripts is chosen as
    the representative sample.
    """
    txs_by_gene = collections.defaultdict(list)
    with open(in_file) as in_handle:
        for rec in SeqIO.parse(in_handle, "fasta"):
            protein_id = rec.id
            header_gene = [p for p in rec.description.split()
                           if p.startswith("gene:")]
            assert len(header_gene) == 1, \
                   "Ensembl gene name not found in header: %s" % rec.description
            (_, gene) = header_gene[0].split(":")
            txs_by_gene[gene].append((len(rec.seq), protein_id))
    final_list = []
    for choices in txs_by_gene.values():
        choices.sort(reverse=True)
        final_list.append(choices[0][1])
    print final_list[:10]
    return set(final_list)

def files_to_filter(base_name):
    base_name = base_name.replace(" ", "_")
    exts = ["ids.tsv", "scores.tsv"]
    fnames_find = ["%s-%s" % (base_name, e) for e in exts]
    fnames = [f for f in fnames_find if os.path.exists(f)]
    if len(fnames) == 0:
        raise ValueError("Did not find files to filter: %s" % fnames_find)
    return fnames

if __name__ == "__main__":
    parser = OptionParser()
    options, args = parser.parse_args()
    if len(args) != 1:
        print __doc__
        sys.exit()
    main(args[0])
