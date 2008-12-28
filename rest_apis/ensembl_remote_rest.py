"""Provide a remote REST-like API interface to Ensembl.

This provides a wrapping around screen-scraping of Ensembl
with the goal of retrieving comparative genomics information.

Usage:
    ensembl_remote_rest.py Organism Gene_id
Example:
    ensembl_remote_rest.py Homo_sapiens ENSG00000173894
"""
from __future__ import with_statement
import sys
import re
import urllib2
import os
import time

from BeautifulSoup import BeautifulSoup
from Bio import SeqIO
import newick
import networkx
    
#gene_ids = [("Homo_sapiens", "ENSG00000173894")]
#gene_ids = [("Homo_sapiens", "ENSG00000168283")]
#gene_ids = [("Homo_sapiens", "ENSG00000183741")]

def main(organism, gene_id):
    write_fasta = False
    cache_dir = os.path.join(os.getcwd(), "cache")
    ensembl_rest = EnsemblComparaRest(cache_dir)
    orthologs = ensembl_rest.orthologs(organism, gene_id)
    compara_tree = ensembl_rest.compara_tree(organism, gene_id)
    compara_tree = '(' + compara_tree[:-1] + ');'
    tree_rec = newick.parse_tree(compara_tree.strip())
    d_vis = DistanceVisitor()
    tree_rec.dfs_traverse(d_vis)
    tree_proteins = [l.identifier for l in tree_rec.leaves]
    orthologs = [(organism, gene_id)] + orthologs
    out_recs = []
    root_id = None
    all_items = []
    for o_organism, o_id in orthologs:
        transcripts = ensembl_rest.transcripts(o_organism, o_id)
        tx, p = [(tx, p) for (tx, p) in transcripts if p in
                tree_proteins][0]
        cur_item = EnsemblComparaTranscript(o_organism, o_id, tx, p)
        if root_id is None:
            root_id = p
        cur_item.distance = networkx.dijkstra_path_length(d_vis.graph,
                "'%s'" % root_id, "'%s'" % p)
        #print o_organism, o_id, p
        cur_item.domains = ensembl_rest.protein_domains(o_organism, o_id,
                tx)
        cur_item.statistics = ensembl_rest.protein_stats(o_organism, o_id,
                tx)
        all_items.append(cur_item)
        if write_fasta:
            out_rec = ensembl_rest.protein_fasta(o_organism, o_id,
                    tx)
            out_rec.id = o_id
            out_rec.description = o_organism
            out_recs.append(out_rec)
    if len(out_recs) > 0:
        with open("%s_%s_orthologs.txt" % (organism, gene_id), "w") as \
                out_handle:
            SeqIO.write(out_recs, out_handle, "fasta")
    analyze_comparative_set(all_items)

def analyze_comparative_set(all_items):
    def distance_cmp(one, two):
        return cmp(one.distance, two.distance)
    all_items.sort(distance_cmp)
    for item in all_items:
        print item.organism, item.distance, item.domains, \
                item.statistics.get('Charge', '').strip()

class EnsemblComparaTranscript:
    """Hold comparative information retrieved from Ensembl on a transcript.
    """
    def __init__(self, organism, g_id, t_id, p_id):
        self.organism = organism
        self.g_id = g_id
        self.t_id = t_id
        self.p_id = p_id
        self.distance = None
        self.domains = []
        self.statistics = {}

class DistanceVisitor(newick.tree.TreeVisitor):
    def __init__(self):
        self.graph = networkx.Graph()
        
    def pre_visit_edge(self, src, b, l, dest):
        self.graph.add_edge(repr(src), repr(dest), l)

class EnsemblComparaRest:
    """Provide a REST-like API interface to Ensembl.
    """
    def __init__(self, cache_dir):
        self._base_url = "http://www.ensembl.org"
        self._cache_dir = cache_dir
        if not(os.path.exists(cache_dir)):
            os.makedirs(cache_dir)

    def protein_stats(self, organism, gene_id, tx_id):
        """Retrieve dictionary of statistics for a gene transcript.
        """
        stats = {}
        with self._get_open_handle("Transcript", "ProteinSummary",
                organism, gene_id, tx_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            stats_possibilities = soup.findAll("dl", "summary")
            for stats_check in stats_possibilities:
                stats_assert = stats_check.find("dt", text="Statistics")
                if stats_assert:
                    stats_line = stats_check.find("p")
                    for stats_part in stats_line:
                        if stats_part.find(":") > 0:
                            key, value = stats_part.split(":")
                            stats[key] = value
        return stats

    def protein_domains(self, organism, gene_id, tx_id):
        """Retrieve characterized domains in a gene transcript.
        """
        domains = []
        with self._get_open_handle("Transcript", "Domains", organism,
                gene_id, tx_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            domain_table = soup.find("table", "ss autocenter")
            if domain_table is not None:
                domain_links = domain_table.findAll("a", href =
                        re.compile("interpro"))
                for domain_link in domain_links:
                    domains.append(domain_link.string)
        domains = list(set(domains))
        return domains

    def protein_fasta(self, organism, gene_id, tx_id):
        """Retrieve the fasta sequence for a given gene and transcript.
        """
        final_url = "%s/%s/Transcript/Export?db=core;g=%s;output=fasta;t=%s;"\
                "st=peptide;_format=Text" % (self._base_url, organism,
                        gene_id, tx_id)
        handle = self._safe_open(final_url)
        rec = SeqIO.read(handle, "fasta")
        handle.close()
        return rec

    def transcripts(self, organism, gene_id):
        """Retrieve a list of (transcript, protein) ids for the given gene_id.
        """
        txs = []
        ps = []
        valid_gene_starts = ["EN", "FB", "AA", "AG"]
        with self._get_open_handle("Gene", "Summary", organism,
                gene_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            tx_info = soup.find("table", {"id" : "transcripts"})
            if tx_info is None:
                tx_info = soup.find(True, {"id" : "transcripts_text"})
            #print tx_info
            tx_links = tx_info.findAll("a", 
                    href = re.compile("Transcript/Summary"))
            for tx_link in tx_links:
                if tx_link.string and tx_link.string[:2] in valid_gene_starts:
                    txs.append(tx_link.string)
            p_links = tx_info.findAll("a", 
                    href = re.compile("Transcript/ProteinSummary"))
            for p_link in p_links:
                if p_link.string:
                    ps.append(p_link.string)
        assert len(txs) == len(ps), (organism, gene_id, txs, ps)
        return zip(txs, ps)

    def orthologs(self, organism, gene_id):
        """Retrieve a list of orthologs for the given gene ID.
        """
        orthologs = []
        with self._get_open_handle("Gene", "Compara_Ortholog",
                organism, gene_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            orth_table = soup.find("table", "orthologues")
            orth_links = orth_table.findAll("a", 
                    href = re.compile("Gene/Summary"))
            for orth_link in orth_links:
                href_parts = [x for x in orth_link['href'].split('/') if x]
                orthologs.append((href_parts[0], orth_link.string))
        return orthologs

    def compara_tree(self, organism, gene_id):
        """Retrieve the comparative tree calculated by compara.
        """
        with self._get_open_handle("Component/Gene/Web/ComparaTree", "text",
                organism, gene_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            tree_details = soup.find("pre")
            return tree_details.string

    def _get_open_handle(self, item_type, action, organism, gene_id,
                         tx_id = None):
        full_url = "%s/%s/%s/%s?g=%s" % (self._base_url, organism,
                item_type, action, gene_id)
        if tx_id:
            full_url += ";t=%s" % (tx_id)
        url_parts = [p for p in full_url.split("/") if p]
        cache_file = os.path.join(self._cache_dir, "_".join(url_parts[1:]))
        if not os.path.exists(cache_file):
            #print full_url, cache_file
            in_handle = self._safe_open(full_url)
            with open(cache_file, 'w') as out_handle:
                out_handle.write(in_handle.read())
            in_handle.close()
        return open(cache_file, 'r')

    def _safe_open(self, url):
        while 1:
            try:
                in_handle = urllib2.urlopen(url)
                return in_handle
            except urllib2.URLError, msg:
                print msg
                time.sleep(5)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
