#!/usr/bin/env python
"""Convert Glimmer gene predictions into protein sequences.

http://www.cbcb.umd.edu/software/glimmer/

Usage:
    glimmer_to_proteins.py <glimmer output> <ref fasta>
"""
from __future__ import with_statement
import sys
import os
import operator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(glimmer_file, ref_file):
    with open(ref_file) as in_handle:
        ref_rec = SeqIO.read(in_handle, "fasta")

    base, ext = os.path.splitext(glimmer_file)
    out_file = "%s-proteins.fa" % base
    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(glimmer_file, ref_rec), out_handle, "fasta")

def protein_recs(glimmer_file, ref_rec):
    """Generate protein records
    """
    with open(glimmer_file) as in_handle:
        for gene_num, exons, strand in glimmer_predictions(in_handle):
            seq_exons = []
            for start, end in exons:
                seq_exons.append(ref_rec.seq[start:end])
            gene_seq = reduce(operator.add, seq_exons)
            if strand == '-':
                gene_seq = gene_seq.reverse_complement()
            protein_seq = gene_seq.translate()
            yield SeqRecord(protein_seq, gene_num, "", "")

def glimmer_predictions(in_handle):
    """Parse Glimmer output, generating a exons and strand for each prediction.
    """
    # read the header
    while 1:
        line = in_handle.readline()
        if line.startswith("   #    #"):
            break
    in_handle.readline()
    # read gene predictions one at a time
    cur_exons, cur_gene_num, cur_strand = ([], None, None)
    while 1:
        line = in_handle.readline()
        if not line:
            break
        parts = line.strip().split()
        # new exon
        if len(parts) == 0:
            yield cur_gene_num, cur_exons, cur_strand
            cur_exons, cur_gene_num, cur_strand = ([], None, None)
        else:
            this_gene_num = parts[0]
            this_strand = parts[2]
            this_start = int(parts[4]) - 1 # 1 based
            this_end = int(parts[5])
            if cur_gene_num is None:
                cur_gene_num = this_gene_num
                cur_strand = this_strand
            else:
                assert cur_gene_num == this_gene_num
                assert cur_strand == this_strand
            cur_exons.append((this_start, this_end))
    if len(cur_exons) > 0:
        yield cur_gene_num, cur_exons, cur_strand

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print __doc__
        sys.exit()
    main(*sys.argv[1:])
