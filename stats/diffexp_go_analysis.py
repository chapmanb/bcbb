#!/usr/bin/env python
"""Provide topGO analysis of overrepresented GO annotation terms in a dataset.

Usage:
    stats_go_analysis.py <input CVS> <gene to GO file>
"""
from __future__ import with_statement
import sys
import csv
import collections

import rpy2.robjects as robjects

def main(input_csv, gene_to_go_file):
    gene_pval = 1e-2
    go_pval = 0.2
    go_term_type = "MF"
    topgo_method = "classic" # choice of classic, elim, weight

    with open(input_csv) as in_handle:
        genes_w_pvals = parse_input_csv(in_handle)
    with open(gene_to_go_file) as in_handle:
        gene_to_go, go_to_gene = parse_go_map_file(in_handle, genes_w_pvals)
    if len(gene_to_go) == 0:
        raise ValueError("No GO terms match to input genes. "
              "Check that the identifiers between the input and GO file match.")
    go_terms = run_topGO(genes_w_pvals, gene_to_go, go_term_type,
            gene_pval, go_pval, topgo_method)
    print_go_info(go_terms, go_term_type, go_to_gene)

def print_go_info(go_terms, go_term_type, go_to_gene):
    for final_pval, go_id, go_term in go_terms:
        genes = []
        for check_go in [go_id] + get_go_children(go_id, go_term_type):
            genes.extend(go_to_gene.get(check_go, []))
        genes = sorted(list(set(genes)))
        print "-> %s (%s) : %0.4f" % (go_id, go_term, final_pval)
        for g in genes:
            print g

def get_go_children(go_term, go_term_type):
    """Retrieve all more specific GO children from a starting GO term.
    """
    robjects.r('''
        library(GO.db)
    ''')
    child_map = robjects.r["GO%sCHILDREN" % (go_term_type)]
    children = []
    to_check = [go_term]
    while len(to_check) > 0:
        new_children = []
        for check_term in to_check:
            new_children.extend(list(robjects.r.get(check_term, child_map)))
        new_children = list(set([c for c in new_children if c]))
        children.extend(new_children)
        to_check = new_children
    children = list(set(children))
    return children

def _dict_to_namedvector(init_dict):
    """Call R to create a named vector from an input dictionary.
    """
    return robjects.r.c(**init_dict)

def run_topGO(gene_vals, gene_to_go, go_term_type, gene_pval, go_pval,
        topgo_method):
    """Run topGO, returning a list of pvalues and terms of interest.
    """
    # run topGO with our GO and gene information
    robjects.r('''
        library(topGO)
    ''')
    robjects.r('''
        topDiffGenes = function(allScore) {
          return (allScore < %s)
        }
    ''' % gene_pval)
    params = {"ontology" : go_term_type,
              "annot" : robjects.r["annFUN.gene2GO"],
              "geneSelectionFun" : robjects.r["topDiffGenes"],
              "allGenes" : _dict_to_namedvector(gene_vals),
              "gene2GO" : _dict_to_namedvector(gene_to_go)
              }
    go_data = robjects.r.new("topGOdata", **params)
    results = robjects.r.runTest(go_data, algorithm=topgo_method,
            statistic="fisher")
    scores = robjects.r.score(results)
    num_summarize = min(100, len(scores.names))
    # extract term names from the topGO summary dataframe
    results_table = robjects.r.GenTable(go_data, elimFisher=results,
            orderBy="elimFisher", topNodes=num_summarize)
    print results_table
    GO_ID_INDEX = 0
    TERM_INDEX = 1
    ids_to_terms = dict()
    for index, go_id in enumerate(results_table[GO_ID_INDEX]):
        ids_to_terms[go_id] = results_table[TERM_INDEX][index]
    go_terms = []
    # convert the scores and results information info terms to return
    for index, item in enumerate(scores):
        if item < go_pval:
            go_id = scores.names[index]
            go_terms.append((item, go_id, ids_to_terms.get(go_id, "")))
    go_terms.sort()
    return go_terms

def parse_go_map_file(in_handle, genes_w_pvals):
    gene_list = genes_w_pvals.keys()
    gene_to_go = collections.defaultdict(list)
    go_to_gene = collections.defaultdict(list)
    for line in in_handle:
        parts = line.split("\t")
        gene_id = parts[0]
        go_id = parts[1].strip()
        if gene_id in gene_list:
            gene_to_go[gene_id].append(go_id)
            go_to_gene[go_id].append(gene_id)
    return dict(gene_to_go), dict(go_to_gene)

def parse_input_csv(in_handle):
    reader = csv.reader(in_handle)
    reader.next() # header
    all_genes = dict()
    for (gene_name, _, _, pval) in reader:
        all_genes[gene_name] = float(pval)
    return all_genes

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print __doc__
        sys.exit()
    main(sys.argv[1], sys.argv[2])
