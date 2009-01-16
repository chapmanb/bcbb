#!/usr/bin/env python
"""Classify a set of proteins from an InterPro query based on descriptions.

This uses a download of InterPro IDs and groups them based on metadata retrieved
through Zemanta analysis of the functional descriptions. This will be useful for
well-characterized organisms that contain hand-written descriptions.

Usage:
    interpro_query_classify.py <input tab file> <API key>
"""
from __future__ import with_statement
import sys
import os
import urllib, urllib2
import xml.etree.ElementTree as ET
import time
import simplejson
import shelve
import operator
import collections

import numpy
from Bio import Cluster

def main(target_id, in_file, api_key):
    cache_dir = os.path.join(os.getcwd(), "cache")
    uniprot_retriever = UniprotRestRetrieval(cache_dir)
    cur_db = shelve.open("%s.db" % os.path.splitext(in_file)[0])
    # load the database
    with open(in_file) as in_handle:
        in_handle.readline() # header
        for index, line in enumerate(in_handle):
            uniprot_id = line.split()[0].strip()
            if uniprot_id not in cur_db.keys():
                cur_terms = get_description_terms(uniprot_retriever,
                        uniprot_id, api_key)
                if len(cur_terms) > 0:
                    cur_db[uniprot_id] = cur_terms
    # cluster and print out cluster details
    term_matrix, uniprot_ids = organize_term_array(cur_db)
    cluster_ids, error, nfound = Cluster.kcluster(term_matrix,
            nclusters=10, npass=20, method='a', dist='e')
    cluster_dict = collections.defaultdict(lambda: [])
    for i, cluster_id in enumerate(cluster_ids):
        cluster_dict[cluster_id].append(uniprot_ids[i])
    for cluster_group in cluster_dict.values():
        if target_id in cluster_group:
            for item in cluster_group:
                print item, cur_db[item]
    cur_db.close()
    
def organize_term_array(cur_db):
    """Organize a set of terms into a binary matrix for classification.

    The rows in the final matrix are the database ids, while the columns are
    terms. Each value is 1 if the term is relevant to that ID and 0 otherwise.
    """
    # flatten all terms and get a unique set
    all_terms = reduce(operator.add, cur_db.values())
    term_counts = collections.defaultdict(lambda: 0)
    for term in all_terms:
        term_counts[term] += 1
    all_terms = list(set(all_terms))
    term_matrix = []
    all_ids = []
    for uniprot_id, cur_terms in cur_db.items():
        cur_row = [(1 if t in cur_terms else 0) for t in all_terms]
        term_matrix.append(cur_row)
        all_ids.append(uniprot_id)
    return numpy.array(term_matrix), all_ids

def get_description_terms(retriever, cur_id, api_key):
    metadata = retriever.get_xml_metadata(cur_id)
    if metadata.has_key("function_descr"):
        #print metadata["function_descr"]
        keywords = zemanta_link_kws(metadata["function_descr"], api_key)
        if len(keywords) > 0:
            return keywords
    return []

def zemanta_link_kws(search_text, api_key):
    """Query Zemanta for keywords linked out to wikipedia or freebase.
    """
    gateway = 'http://api.zemanta.com/services/rest/0.0/'
    args = {'method': 'zemanta.suggest',
            'api_key': api_key,
            'text': search_text,
            'return_categories': 'dmoz',
            'return_images': 0,
            'return_rdf_links' : 1,
            'format': 'json'}
    args_enc = urllib.urlencode(args)
    raw_output = urllib2.urlopen(gateway, args_enc).read()
    output = simplejson.loads(raw_output)

    link_kws = []
    for link in output['markup']['links']:
        for target in link['target']:
            if target['type'] in ['wikipedia', 'rdf']:
                link_kws.append(target['title'])
    return list(set(link_kws))

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
                if str(msg).find("404: Not Found") >= 0:
                    self._add_not_found(url)
                    return None
                print msg
                time.sleep(5)

    def _add_not_found(self, url):
        with open(self._not_found_file, 'a') as out_handle:
            out_handle.write("%s\n" % url)
        self._not_found.append(url)

class UniprotRestRetrieval(_BaseCachingRetrieval):
    """Retrieve RDF data from UniProt for proteins of interest.
    """
    def __init__(self, cache_dir):
        _BaseCachingRetrieval.__init__(self, cache_dir)
        self._server = "http://www.uniprot.org"
        self._xml_ns = "{http://uniprot.org/uniprot}"

    def get_xml_metadata(self, uniprot_id):
        """Retrieve data from the UniProt XML for a record.

        XXX This retrieves only a subset of metadata right now. Needs to
        be complete.
        """
        url_base = "%s/uniprot/%s.xml"
        full_url = url_base % (self._server, uniprot_id)
        # check for empty files -- which have been deleted
        with self._get_open_handle(full_url) as in_handle:
            if in_handle.readline() == "":
                return {}
        metadata = {}
        with self._get_open_handle(full_url) as in_handle:
            root = ET.parse(in_handle).getroot()
            metadata = self._get_org_metadata(root, metadata)
            metadata = self._get_interpro_metadata(root, metadata)
            metadata = self._get_function_metadata(root, metadata)
        return metadata

    def _get_org_metadata(self, root, metadata):
        """Retrieve the organism information from UniProt XML.
        """
        org = root.find("%sentry/%sorganism" % (self._xml_ns, self._xml_ns))
        for org_node in org:
            if org_node.tag == "%sname" % self._xml_ns:
                if org_node.attrib["type"] == "scientific":
                    metadata["org_scientific_name"] = org_node.text
                elif org_node.attrib["type"] == "common":
                    metadata["org_common_name"] = org_node.text
            elif org_node.tag == "%slineage" % self._xml_ns:
                metadata["org_lineage"] = [n.text for n in org_node]
        return metadata

    def _get_interpro_metadata(self, root, metadata):
        """Retrieve InterPro domains present in the protein.
        """
        db_refs = root.findall("%sentry/%sdbReference" % (self._xml_ns,
            self._xml_ns))
        all_refs = []
        for db_ref in db_refs:
            if db_ref.attrib["type"] in ["InterPro"]:
                all_refs.append("%s:%s" % (db_ref.attrib["type"],
                    db_ref.attrib["id"]))
        if len(all_refs) > 0:
            metadata["db_refs"] = all_refs
        return metadata

    def _get_function_metadata(self, root, metadata):
        """Retrieve an InterPro function description.
        """
        comments = root.findall("%sentry/%scomment" % (self._xml_ns,
            self._xml_ns))
        for comment in comments:
            if comment.attrib["type"] in ["function"]:
                for comment_node in comment:
                    if comment_node.tag == "%stext" % (self._xml_ns):
                        metadata["function_descr"] = comment_node.text
        return metadata

    def get_rdf_metadata(self, uniprot_id):
        """Retrieve RDF metadata for the given UniProt accession.

        XXX Not finished. XML parsing looks to be more straightforward
        """
        from rdflib import ConjunctiveGraph as Graph
        url_base = "%s/uniprot/%s.rdf"
        full_url = url_base % (self._server, uniprot_id)
        graph = Graph()
        with self._get_open_handle(full_url) as in_handle:
            graph.parse(in_handle)
        main_subject = [s for s in list(set(graph.subjects())) if
                s.split('/')[-1] == uniprot_id][0]
        for sub, pred, obj in graph:
            print sub, pred, obj

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1], sys.argv[2], sys.argv[3])
