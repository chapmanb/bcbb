#!/usr/bin/env python
"""Characterize proteins available in InterPro containing a specific domain.

This script takes as input an InterPro domain of interest, pulls down
proteins containing that domain, and provides an output table of those
proteins plus characteristics for downstream analysis.

Usage:
    interpro_domain_summary.py <IPRnumber>
"""
from __future__ import with_statement
import sys
import os
import urllib2
import subprocess
import xml.etree.ElementTree as ET
import shelve
import time

from Bio import SeqIO
from Bio.Emboss import Applications

def main(ipr_number):
    cache_dir = os.path.join(os.getcwd(), "cache")
    db_dir = os.path.join(os.getcwd(), "db")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    interpro_retriever = InterproRestRetrieval(cache_dir)
    uniprot_retriever = UniprotRestRetrieval(cache_dir)
    charge_calc = ProteinChargeCalculator(cache_dir)
    cur_db = shelve.open(os.path.join(db_dir, ipr_number))
    seq_recs = interpro_retriever.get_interpro_records(ipr_number)
    for seq_rec in seq_recs:
        uniprot_id = seq_rec.id.split("|")[0]
        metadata = uniprot_retriever.get_xml_metadata(uniprot_id)
        if (metadata.has_key("org_lineage") and 
                "Metazoa" in metadata["org_lineage"]):
            metadata["seq"] = seq_rec.seq.data
            metadata["charge"] = charge_calc.get_neutral_charge(seq_rec)
            if metadata.has_key("function_descr"):
                print uniprot_id, metadata
            cur_db[uniprot_id] = metadata
    cur_db.close()

class ProteinChargeCalculator:
    """Calculate protein charge using the emboss iep program.
    """
    def __init__(self, tmp_dir):
        self._tmp_dir = tmp_dir

    def get_neutral_charge(self, aa_rec):
        """Get the charge of the provided protein at pH 7.
        """
        in_file = os.path.join(self._tmp_dir, "iep_in.txt")
        out_file = os.path.join(self._tmp_dir, "iep_out.txt")
        try:
            charge = self._calc_neutral_charge(aa_rec, in_file, out_file,
                    '7.00')
        finally:
            for fname in [in_file, out_file]:
                if os.path.exists(fname):
                    os.remove(fname)
        return charge

    def _calc_neutral_charge(self, aa_rec, in_file, out_file, ph_target):
        with open(in_file, 'w') as in_handle:
            SeqIO.write([aa_rec], in_handle, "fasta")
        iep_cl = Applications.IepCommandline()
        iep_cl.set_parameter('-sequence', in_file)
        iep_cl.set_parameter('-outfile', out_file)
        child = subprocess.Popen(str(iep_cl).split(), stdout =
                subprocess.PIPE, stderr = subprocess.PIPE)
        child.wait()
        return self._parse_out_file(out_file, ph_target)

    def _parse_out_file(self, out_file, ph_target):
        with open(out_file) as out_handle:
            for line in out_handle:
                parts = [p for p in line.strip().split() if p]
                if len(parts) == 3 and parts[0] == ph_target:
                    return float(parts[2])
        raise ValueError('Did not find a charge')

class _BaseCachingRetrieval:
    """Provide a base class for web retrieval with local file caching.
    """
    def __init__(self, cache_dir):
        self._cache_dir = cache_dir
        if not(os.path.exists(cache_dir)):
            os.makedirs(cache_dir)

    def _get_open_handle(self, full_url):
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

class InterproRestRetrieval(_BaseCachingRetrieval):
    """Interface to retrieving data via REST interfaces at EBI interpro.
    """
    def __init__(self, cache_dir):
        _BaseCachingRetrieval.__init__(self, cache_dir)
        self._server = "http://www.ebi.ac.uk"

    def get_interpro_records(self, ipr_number):
        """Retrieve protein Biopython SeqRecords with the given domain.
        """
        url_base = "%s/interpro/ISpy?ipr=%s&mode=fasta"
        full_url = url_base % (self._server, ipr_number)
        recs = []
        with self._get_open_handle(full_url) as in_handle:
            for rec in SeqIO.parse(in_handle, "fasta"):
                recs.append(rec)
        return recs

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1])
