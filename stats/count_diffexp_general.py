#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.

http://bioconductor.org/packages/2.5/bioc/html/edgeR.html

Usage:
    count_diffexp.py <count_file>
"""
import os
import sys
import csv
import collections

import numpy
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

def main(count_file):
    base, ext = os.path.splitext(count_file)
    outfile = "%s-diffs.csv" % (base)
    counts, all_regions, conditions, groups, sizes = read_count_file(count_file)
    data, regions, sizes = edger_matrix(counts, conditions, all_regions)
    probs = run_edger(data, groups, sizes, regions)
    write_outfile(outfile, regions, conditions, counts, probs, sizes)

def write_outfile(outfile, genes, conditions, work_counts, probs, sizes):
    with open(outfile, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["Region"] +
                ["%s count" % c for c in conditions] + ["edgeR p-value"])
        writer.writerow(["total"] + [str(s) for s in sizes])
        out_info = []
        for i, gene in enumerate(genes):
            counts = [int(work_counts[c][gene]) for c in conditions]
            out_info.append((probs[i], [gene] + counts))
        out_info.sort()
        [writer.writerow(start + [prob]) for prob, start in out_info]

def run_edger(data, groups, sizes, regions):
    """Call edgeR in R and organize the resulting differential expressed genes.
    """
    robjects.r('''
        library(edgeR)
    ''')
    # find the version we are running -- check for edgeR exactTest function
    try:
        robjects.r["exactTest"]
    except LookupError:
        raise ValueError("Need edgeR 1.3+ to run analysis.")
    params = {'group' : numpy.array(groups), 'lib.size' : sizes}
    dgelist = robjects.r.DGEList(data, **params)
    # perform Poisson adjustment and assignment as recommended in the manual
    robjects.globalEnv['dP'] = dgelist
    # if we have replicates, can estimate common and tagwise dispersion
    if len(groups) > 2:
        robjects.r('''
          dP <- estimateCommonDisp(dP)
          prior.weight <- estimateSmoothing(dP)
          dP <- estimateTagwiseDisp(dP, prior.n=10)
        ''')
    # otherwise use a Poisson distribution estimation (Section 9 of manual)
    else:
        robjects.r('''
            msP <- de4DGE(dP, doPoisson = TRUE)
            dP$pseudo.alt <- msP$pseudo
            dP$common.dispersion <- 1e-06
            dP$conc <- msP$conc
            dP$common.lib.size <- msP$M
        ''')
    dgelist = robjects.globalEnv['dP']
    de = robjects.r.exactTest(dgelist)
    tag_table = robjects.r.topTags(de, n=len(regions))[0]
    print robjects.r.head(tag_table)
    indexes = [int(t.replace("tag.", "")) - 1 for t in tag_table.rownames()]
    # can retrieve either raw or adjusted p-values
    #pvals = list(tags.r['p.value'][0])
    pvals = list(tag_table.r['PValue'][0])

    assert len(indexes) == len(pvals)
    pvals_w_index = zip(indexes, pvals)
    pvals_w_index.sort()
    assert len(pvals_w_index) == len(indexes)
    return [p for i,p in pvals_w_index]

def edger_matrix(work_counts, conditions, regions):
    """Count matrices for input into edgeR differential expression analysis.
    """
    data = []
    final_regions = []
    for r in regions:
        cur_row = [int(work_counts[c][r]) for c in conditions]
        if sum(cur_row) > 0:
            final_regions.append(r)
            data.append(cur_row)
    sizes = numpy.sum(data, axis=0)
    return numpy.array(data), final_regions, sizes

def read_count_file(in_file):
    """Read count information from a simple CSV file into a dictionary.
    """
    counts = collections.defaultdict(dict)
    regions = []
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        conditions = reader.next()[1:]
        groups = [int(p) for p in reader.next()[1:]]
        totals = [int(p) for p in reader.next()[1:]]
        for parts in reader:
            region_name = parts[0]
            regions.append(region_name)
            region_counts = [float(x) for x in parts[1:]]
            for ci, condition in enumerate(conditions):
                counts[condition][region_name] = region_counts[ci]
    return dict(counts), regions, conditions, groups, totals

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print __doc__
        sys.exit()
    main(sys.argv[1])
