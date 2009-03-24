#!/usr/bin/env python
"""Characterize proteins available in InterPro containing a specific domain.

This script takes as input an InterPro domain of interest, pulls down
proteins containing that domain, and provides an output table of those
proteins plus characteristics for downstream analysis.

Usage:
    interpro_domain_summary.py <IPRnumber> <domain description>
"""
from __future__ import with_statement
import sys
import os
import urllib2
import subprocess
import xml.etree.ElementTree as ET
import shelve
import time
import collections

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss import Applications
from Bio.SeqUtils import IsoelectricPoint
from Bio.SeqUtils import ProtParam

org_includes = [
        'Homo sapiens',
        'Mus musculus',
        'Rattus norvegicus',
        'Bos taurus',
        'Gallus gallus',
        'Xenopus laevis',
        'Danio rerio',
        'Tetraodon nigroviridis',
        'Ciona intestinalis',
        'Branchiostoma floridae',
        'Caenorhabditis elegans',
        'Anopheles gambiae',
        'Drosophila melanogaster',
        'Nematostella vectensis'
]

extra_info = ['Ciona savignyi', 'Ciona intestinalis',
        'Strongylocentrotus purpuratus']

def main(ipr_number, domain_description):
    aa_window = 75
    charge_thresh = 10.2
    cache_dir = os.path.join(os.getcwd(), "cache")
    db_dir = os.path.join(os.getcwd(), "db")
    domain_dir = os.path.join(os.getcwd(), "pfam_domains")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    interpro_retriever = InterproRestRetrieval(cache_dir)
    uniprot_retriever = UniprotRestRetrieval(cache_dir)
    uniref_retriever = UniRefRetrieval(cache_dir)
    string_retriever = StringRetrieval(cache_dir)
    charge_calc = ProteinChargeCalculator(cache_dir)
    hmmsearch_parser = SimpleHmmsearchDomainParser(std_name_parser)
    cur_db = shelve.open(os.path.join(db_dir, ipr_number))
    seq_recs = interpro_retriever.get_interpro_records(ipr_number)
    with open(os.path.join(domain_dir, 
              "%s-domains.hmmsearch" % ipr_number)) as in_handle:
        hmm_locations = hmmsearch_parser.domain_locations(in_handle)
    # retrieve filtering information based on HMM searches
    with open(os.path.join(domain_dir, 
              "%s-filter1.hmmsearch" % ipr_number)) as in_handle:
        filter_hmm = hmmsearch_parser.domain_locations(in_handle)
    uniref_data = {}
    all_children = []
    all_recs = []
    for seq_rec in seq_recs:
        uniprot_id = std_name_parser(seq_rec.id)
        metadata = uniprot_retriever.get_xml_metadata(uniprot_id,
                domain_description)
        if (metadata.has_key("org_lineage") and 
                "Metazoa" in metadata["org_lineage"] and
                (len(org_includes) == 0 or 
                 metadata["org_scientific_name"] in org_includes)):
            all_recs.append(seq_rec)
            if metadata["org_scientific_name"] in extra_info:
                print metadata["org_scientific_name"], uniprot_id, \
                        passes_filters(uniprot_id, filter_hmm)
            if passes_filters(uniprot_id, filter_hmm):
                if not metadata.has_key("domain_positions"):
                    positions = []
                    for start, end in hmm_locations.get(uniprot_id, []):
                        positions.append(start)
                        positions.append(end)
                    metadata["domain_positions"] = positions
                if len(metadata["domain_positions"]) > 0:
                    interactors = string_retriever.get_string_interactions(
                            uniprot_id)
                    if len(interactors) > 0:
                        metadata["string_interactors"] = interactors
                    uniref_info = uniref_retriever.get_90_group(uniprot_id)
                    for org, vals in uniref_info.items():
                        uniref_data[vals[0]] = vals[1:]
                        all_children.extend(vals[1:])
                    metadata["seq"] = seq_rec.seq.data
                    metadata["charge"] = charge_calc.get_neutral_charge(seq_rec)
                    metadata["charge_region"] = \
                            charge_calc.get_region_charge_percent(
                            seq_rec, aa_window, charge_thresh)
                    #if metadata.has_key("function_descr"):
                    #    print uniprot_id, metadata
                    cur_db[uniprot_id] = metadata
    # add information about whether an item is a child or parent in
    # a uniref group based on sequence similarity
    for uniprot_id in cur_db.keys():
        metadata = cur_db[uniprot_id]
        try:
            uniref_children = uniref_data[uniprot_id]
            metadata["uniref_children"] = uniref_children
        except KeyError:
            if uniprot_id in all_children:
                metadata["is_uniref_child"] = "yes"
        cur_db[uniprot_id] = metadata
    with open(os.path.join(domain_dir, "%s.fa" % ipr_number),
            "w") as out_handle:
        SeqIO.write(all_recs, out_handle, "fasta")
    cur_db.close()

def passes_filters(uniprot_id, filter_hmm):
    """Determine if a uniprot_id contains domains based on HMM search results.

    This is meant as a general filter to retrieve only items contained defined
    domains regions.
    """
    return len(filter_hmm.get(uniprot_id, [])) > 0

def std_name_parser(name):
    return name.split('|')[0]

class SimpleHmmsearchDomainParser:
    """Simple hmmsearch parser which pulls out domain locations as a dictionary.
    """
    def __init__(self, name_parser=None):
        """Initialize with a function to parse out name information.

        The name parser is a function that will parse out the final name 
        from that listed in the Hmm
        """
        # define a do nothing function if one is not supplied
        if name_parser is None:
            name_parser = lambda x: x
        self._name_parser = name_parser

    def domain_locations(self, in_handle):
        """Retrieve a dictionary of domain locations from a Hmmsearch report.
        """
        while 1:
            line = in_handle.readline()
            if line.find('Parsed for domains:') == 0:
                break
        in_handle.readline()
        in_handle.readline()
        domain_locs = collections.defaultdict(lambda: []) 
        while 1:
            line = in_handle.readline().strip()
            if not line:
                break
            name, start, end = self._parse_domain_line(line)
            domain_locs[name].append((start, end))
        return dict(domain_locs)

    def _parse_domain_line(self, line):
        """Retrieve domain information from the Hmmsearch report.

        Q28HD9|Q28HD9_XENTR   1/2  21    70 ..   1    69 []    90.0 7.7e-25
        """
        parts = [x.strip() for x in line.split() if x.strip()]
        name = self._name_parser(parts[0])
        start = int(parts[2]) - 1
        end = int(parts[3])
        return name, start, end

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

    def get_region_charge_percent(self, aa_rec, aa_window, charge_thresh):
        """Retrieve the percent of windows in a sequence above a given threshold
        """
        region_charges = self._calc_region_charges(aa_rec.seq, aa_window)
        above_thresh = [c for c in region_charges if c >= charge_thresh]
        if len(region_charges) > 0:
            return float(len(above_thresh)) / float(len(region_charges))
        else:
            return 0.0

    def _calc_region_charges(self, seq, cur_window):
        """Perform calculation of charges via isoelectric points for a sequence.
        """
        # internal small regions, so do not deal with C and N terminal charges
        IsoelectricPoint.pKcterminal = {}
        IsoelectricPoint.pKnterminal = {}
        cur_pos = 0
        region_charges = []
        while cur_pos < len(seq) - cur_window:
            cur_seq = seq[cur_pos:cur_pos + cur_window]
            prot_analysis = ProtParam.ProteinAnalysis(str(cur_seq))
            ie_calc = IsoelectricPoint.IsoelectricPoint(cur_seq,
                    prot_analysis.count_amino_acids())
            region_charges.append(ie_calc.pi())
            cur_pos += 1
        return region_charges

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

class StringRetrieval(_BaseCachingRetrieval):
    """Retrieve protein interaction data from STRING.
    """
    def __init__(self, cache_dir):
        _BaseCachingRetrieval.__init__(self, cache_dir)
        self._server = "http://string.embl.de"

    def get_string_interactions(self, protein_id):
        """Retrieve interaction data from STRING.
        """
        url_base = "%s/api/tsv/interactors?identifier=%s"
        full_url = url_base % (self._server, protein_id)
        in_handle = self._get_open_handle(full_url)
        # 400 error message
        if in_handle is None:
            return []
        with in_handle:
            lines = in_handle.read().split('\n')
            # if we get an error, no interactions found
            if lines[0].find('ErrorMessage') >= 0:
                return []
            else:
                return [l for l in lines[1:] if l]

class UniRefRetrieval(_BaseCachingRetrieval):
    """Retrieve UniRef clusters for provided
    """
    def __init__(self, cache_dir):
        _BaseCachingRetrieval.__init__(self, cache_dir)
        self._server = "http://www.uniprot.org"
        self._xml_ns = "{http://uniprot.org/uniref}"

    def get_90_group(self, uniprot_id):
        """Retrieve the UniProt 90 cluster related to the given uniprot_id.
        """
        url_base = "%s/uniref/UniRef90_%s.xml"
        full_url = url_base % (self._server, uniprot_id)
        group_ids = collections.defaultdict(lambda: [])
        in_handle = self._get_open_handle(full_url)
        if in_handle:
            with in_handle:
                root = ET.parse(in_handle).getroot()
                db_refs = root.findall(
                            "%sentry/%srepresentativeMember/%sdbReference" % (
                            (self._xml_ns,) * 3)) + \
                          root.findall("%sentry/%smember/%sdbReference" % (
                            (self._xml_ns,) * 3))
                for db_ref in db_refs:
                    if db_ref.attrib["type"] in ["UniProtKB ID"]:
                        cur_id, cur_org = (None, None)
                        for prop in db_ref:
                            if prop.attrib["type"] in ["UniProtKB accession"]:
                                cur_id = prop.attrib["value"]
                            if prop.attrib["type"] in ["source organism"]:
                                cur_org = prop.attrib["value"]
                        if cur_id and cur_org:
                            cur_org = " ".join(cur_org.split()[:2])
                            group_ids[cur_org].append(cur_id)
        return self._filter_groups(dict(group_ids))

    def _filter_groups(self, group_ids):
        """Return only groups that contain more than a single identifier.
        """
        final_group_ids = {}
        for org, vals in group_ids.items():
            if len(vals) > 1:
                final_group_ids[org] = vals
        return final_group_ids

class UniprotRestRetrieval(_BaseCachingRetrieval):
    """Retrieve RDF data from UniProt for proteins of interest.
    """
    def __init__(self, cache_dir):
        _BaseCachingRetrieval.__init__(self, cache_dir)
        self._server = "http://www.uniprot.org"
        self._xml_ns = "{http://uniprot.org/uniprot}"

    def get_xml_metadata(self, uniprot_id, domain_description):
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
            metadata = self._get_genename_metadata(root, metadata)
            metadata = self._get_org_metadata(root, metadata)
            metadata = self._get_dbref_metadata(root, metadata)
            metadata = self._get_function_metadata(root, metadata)
            metadata = self._get_domain_position_metadata(domain_description,
                    root, metadata)
        return metadata

    def _get_genename_metadata(self, root, metadata):
        """Retrieve gene names for this uniprot record.
        """
        gene = root.find("%sentry/%sgene" % (self._xml_ns, self._xml_ns))
        gene_names = []
        if gene:
            for gene_node in gene:
                if gene_node.tag == "%sname" % self._xml_ns:
                    gene_names.append(gene_node.text)
        if len(gene_names) > 0:
            metadata["gene_names"] = gene_names
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

    def _get_dbref_metadata(self, root, metadata):
        """Retrieve InterPro and Ensembl DB xrefs present in the protein.

        XXX This needs to be generalized badly.
        """
        db_refs = root.findall("%sentry/%sdbReference" % (self._xml_ns,
            self._xml_ns))
        all_refs = []
        ensembl_refs = []
        for db_ref in db_refs:
            if db_ref.attrib["type"] in ["InterPro"]:
                all_refs.append("%s:%s" % (db_ref.attrib["type"],
                    db_ref.attrib["id"]))
            elif db_ref.attrib["type"] in ["Ensembl"]:
                ensembl_refs.append(db_ref.attrib["id"])
        if len(all_refs) > 0:
            metadata["db_refs"] = all_refs
        if len(ensembl_refs) > 0:
            metadata["db_refs_ensembl"] = ensembl_refs
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

    def _get_domain_position_metadata(self, domain_descr, root, metadata):
        """Find the positions of our domain of interest in the protein.
        """
        features = root.findall("%sentry/%sfeature" % (self._xml_ns,
            self._xml_ns))
        positions = []
        for feature in features:
            if (feature.attrib["type"] in ["domain", "zinc finger region"] and
                feature.attrib["description"].find(domain_descr) >= 0):
                location = feature[0]
                for loc_info in location:
                    cur_pos = int(loc_info.attrib['position'])
                    # convert to 0-based coordinates
                    if loc_info.tag.find("begin") >= 0:
                        cur_pos -= 1
                    positions.append(cur_pos)
        if len(positions) > 0:
            metadata["domain_positions"] = positions
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
    if len(sys.argv) != 3:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(sys.argv[1], sys.argv[2])
