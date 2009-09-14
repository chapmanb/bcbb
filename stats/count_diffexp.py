#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.

http://bioconductor.org/packages/2.5/bioc/html/edgeR.html

Usage:
    count_diffexp.py <count_file>
"""
import math
import os
import sys
import csv
import collections
import itertools

import numpy
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

def main(count_file):
    base, ext = os.path.splitext(count_file)
    outfile = "%s-diffs.csv" % (base)
    counts = read_count_file(count_file)
    data, groups, sizes, conditions, genes = edger_matrices(counts)
    probs = run_edger(data, groups, sizes, genes)
    write_outfile(outfile, genes, conditions, counts, probs)

def write_outfile(outfile, genes, conditions, work_counts, probs):
    with open(outfile, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["Region"] +
                ["%s count" % c for c in conditions] + ["edgeR p-value"])
        out_info = []
        for i, gene in enumerate(genes):
            counts = [int(work_counts[c][gene]) for c in conditions]
            out_info.append((probs[i], [gene] + counts))
        out_info.sort()
        [writer.writerow(start + [prob]) for prob, start in out_info]

def run_edger(data, groups, sizes, genes):
    """Call edgeR in R and organize the resulting differential expressed genes.
    """
    robjects.r('''
        library(edgeR)
    ''')
    params = {'group' : groups, 'lib.size' : sizes}
    dgelist = robjects.r.DGEList(data, **params)
    ms = robjects.r.deDGE(dgelist, doPoisson=True)
    tags = robjects.r.topTags(ms, pair=groups, n=len(genes))
    indexes = [int(t) - 1 for t in tags.rownames()]
    pvals = list(tags[-2])
    assert len(indexes) == len(pvals)
    pvals_w_index = zip(indexes, pvals)
    pvals_w_index.sort()
    assert len(pvals_w_index) == len(indexes)
    return [p for i,p in pvals_w_index]

def get_conditions_and_genes(work_counts): 
    conditions = work_counts.keys()
    conditions.sort()
    all_genes = []
    for c in conditions:
        all_genes.extend(work_counts[c].keys())
    all_genes = list(set(all_genes))
    all_genes.sort()

    sizes = [work_counts[c]["Total"] for c in conditions]
    all_genes.remove("Total")
    return conditions, all_genes, sizes
    
def edger_matrices(work_counts):
    """Retrieve matrices for input into edgeR differential expression analysis.
    """
    conditions, all_genes, sizes = get_conditions_and_genes(work_counts)
    assert len(sizes) == 2
    groups = [1, 2]
    data = []
    final_genes = []
    for g in all_genes:
        cur_row = [int(work_counts[c][g]) for c in conditions]
        if sum(cur_row) > 0:
            data.append(cur_row)
            final_genes.append(g)
    return (numpy.array(data), numpy.array(groups), numpy.array(sizes),
            conditions, final_genes)

def read_count_file(in_file):
    """Read count information from a simple CSV file into a dictionary.
    """
    counts = collections.defaultdict(dict)
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        header = reader.next()
        conditions = header[1:]
        for parts in reader:
            region_name = parts[0]
            region_counts = [float(x) for x in parts[1:]]
            for ci, condition in enumerate(conditions):
                counts[condition][region_name] = region_counts[ci]
    return dict(counts)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print __doc__
        sys.exit()
    main(sys.argv[1])
