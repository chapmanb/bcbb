#!/usr/bin/env python
"""Cluster together proteins based on genomic location and duplication history.

Given a list of identified proteins of interest, this looks for similarities
based on either genomic location (parsed from Ensembl) or duplication (parsed
from Ensembl compara).

Usage:
    loc_dup_cluster.py <list of uniprot IDs>
"""
from __future__ import with_statement
import sys
import os
import glob
import shelve
import collections
import urllib2
import re
import time

from BeautifulSoup import BeautifulSoup

def main(id_file):
    distance_thresh = 1000000 # locations within 1 megabase
    db_dir = os.path.join(os.getcwd(), "db")
    cache_dir = os.path.join(os.getcwd(), "cache")
    ensembl_retriever = EnsemblComparaRest(cache_dir)
    cur_dbs = get_available_dbs(db_dir)
    all_ids, id_groups = read_in_file(id_file)
    org_locations = collections.defaultdict(lambda: [])
    org_info = collections.defaultdict(lambda: [])
    loc_info = dict()
    dup_info = dict()
    ensembl_to_uniprot = dict()
    for cur_id in all_ids:
        cur_rec = None
        for db in cur_dbs:
            try:
                cur_rec = db[cur_id]
                break
            except KeyError:
                pass
        assert cur_rec is not None, cur_id
        cur_ensembl_org = cur_rec["org_scientific_name"].replace(" ", "_")
        for ensembl_id in cur_rec.get("db_refs_ensembl", []):
            paralogs = ensembl_retriever.paralogs(cur_ensembl_org, ensembl_id)
            chromosome, start, end = ensembl_retriever.location(
                    cur_ensembl_org, ensembl_id)
            org_locations[cur_rec["org_scientific_name"]].append((chromosome,
                start, end))
            dup_info[ensembl_id] = paralogs
            ensembl_to_uniprot[ensembl_id] = cur_id
            org_info[cur_id].append((cur_rec["org_scientific_name"],
                chromosome, start, end))
            loc_info[ensembl_id] = (cur_rec["org_scientific_name"],
                chromosome, start, end)
    dup_groups = examine_paralogs(dup_info, ensembl_to_uniprot)
    print '=> Duplications'
    for dup_group in dup_groups:
        print '---------'
        for uniprot_id in dup_group:
            print uniprot_id, id_groups[uniprot_id], org_info[uniprot_id]
    loc_groups = examine_location(loc_info, ensembl_to_uniprot, distance_thresh)
    print '=> Locations'
    for loc_group in loc_groups:
        print '---------'
        for uniprot_id in loc_group:
            print uniprot_id, id_groups[uniprot_id], org_info[uniprot_id]

def read_in_file(id_file):
    all_ids = []
    id_groups = {}
    with open(id_file) as in_handle:
        for line in in_handle:
            cur_id, id_group = line.strip().split()
            all_ids.append(cur_id)
            id_groups[cur_id] = id_group
    return all_ids, id_groups

def location_distance(loc_one, loc_two):
    """Determine the distance between two locations.

    Locations not in the same organism and not on the same chromosome get a
    huge maxint number to indicate they are not close. Otherwise, this is the
    distance on the chromosome.
    """
    if loc_one[:2] != loc_two[:2] or loc_one == loc_two:
        return sys.maxint
    else:
        return max(abs(loc_one[3] - loc_two[2]),
                   abs(loc_two[3] - loc_one[2]))

def examine_location(loc_info, ensembl_to_uniprot, distance_thresh):
    """Look for gene pairs which are structurally located near each other.
    """
    loc_close = collections.defaultdict(lambda: [])
    for uniprot_id, loc_one in loc_info.items():
        for cmp_id, loc_two in loc_info.items():
            if location_distance(loc_one, loc_two) <= distance_thresh:
                loc_close[uniprot_id].append(cmp_id)
    return examine_paralogs(loc_close, ensembl_to_uniprot)

def examine_paralogs(dup_info, ensembl_to_uniprot):
    """Extract details on paralog groups identified in Ensembl.
    """
    cur_groups = []
    all_base = dup_info.keys()
    for base_id, dup_ids in dup_info.items():
        overlap = set(dup_ids) & set(all_base)
        if len(overlap) > 0:
            new_group = set([ensembl_to_uniprot[x] for x in overlap |
                set([base_id])])
            is_unique = True
            for exist_i, exist_group in enumerate(cur_groups):
                if len(new_group & exist_group) > 0:
                    update_group = new_group | exist_group
                    cur_groups[exist_i] = update_group
                    is_unique = False
                    break
            if is_unique:
                cur_groups.append(new_group)
    return [list(g) for g in cur_groups]

def get_available_dbs(db_dir):
    all_dbs = []
    for fname in glob.glob(os.path.join(db_dir, '*.db')):
        base, ext = os.path.splitext(fname)
        all_dbs.append(shelve.open(base))
    return all_dbs

class _BaseCachingRetrieval:
    """Provide a base class for web retrieval with local file caching.
    """
    def __init__(self, cache_dir):
        self._cache_dir = cache_dir
        if not(os.path.exists(cache_dir)):
            os.makedirs(cache_dir)
        # cache 404 errors so we don't call the page multiple times
        self._not_found_file = os.path.join(self._cache_dir,
                '404_not_found.txt')
        self._not_found = []
        if os.path.exists(self._not_found_file):
            with open(self._not_found_file) as in_handle:
                self._not_found = in_handle.read().split()

    def _get_open_handle(self, full_url):
        if full_url in self._not_found:
            return None
        url_parts = [p for p in full_url.split("/") if p]
        cache_file = os.path.join(self._cache_dir, "_".join(url_parts[1:]))
        if not os.path.exists(cache_file):
            print full_url
            #print full_url, cache_file
            in_handle = self._safe_open(full_url)
            if in_handle is None:
                return None
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
                if (str(msg).find("404: Not Found") >= 0 or
                    str(msg).find("HTTP Error 400: Bad Request") >= 0):
                    self._add_not_found(url)
                    return None
                print msg, url
                time.sleep(5)

    def _add_not_found(self, url):
        with open(self._not_found_file, 'a') as out_handle:
            out_handle.write("%s\n" % url)
        self._not_found.append(url)

class EnsemblComparaRest(_BaseCachingRetrieval):
    """Provide a REST-like API interface to Ensembl.
    """
    def __init__(self, cache_dir):
        self._base_url = "http://www.ensembl.org"
        _BaseCachingRetrieval.__init__(self, cache_dir)

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

    def location(self, organism, gene_id):
        """Retrieve location information for the given organism and gene ID.
        """
        with self._get_open_handle("Gene", "Summary",
                organism, gene_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            loc_tab = soup.find("dd", id="tab_location")
            link = loc_tab.find("a")
            path, attrs = urllib2.splitattr(link["href"])
            for attr in attrs:
                if attr.find("r=") == 0:
                    key, val = attr.split("=")
                    chrom, location = val.split(":")
                    start, end = location.split("-")
                    return chrom, int(start), int(end)
            raise ValueError("Did not find location: %s" % link)

    def paralogs(self, organism, gene_id):
        """Retrieve a list of paralogs for the given gene ID.
        """
        paralogs = []
        with self._get_open_handle("Gene", "Compara_Paralog",
                organism, gene_id) as in_handle:
            soup = BeautifulSoup(in_handle)
            tables = soup.findAll("table")
            for table in tables:
                paralog_links = table.findAll("a",
                        href=re.compile("Gene/Compara_Paralog?"))
                for para_link in paralog_links:
                    if para_link.string not in ["Align"]:
                        paralogs.append(para_link.string)
        return paralogs

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
        return _BaseCachingRetrieval._get_open_handle(self, full_url)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1])
