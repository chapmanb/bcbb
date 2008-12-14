"""Provide a (semi) automated mapping between GenBank and ontologies.

The goal is to provide an ontology namespace and term for each item
of a standard GenBank file.
"""
from __future__ import with_statement
import os
import collections
import csv

from rdflib.Graph import Graph
#from Mgh.SemanticLiterature import sparta
import sparta

def main(ft_file, so_ft_map_file, so_file, dc_file, rdfs_file):
    # -- get all of the keys we need to map from GenBank
    # pulled by hand from the Biopython parser
    header_keys = ['ACCESSION', 'AUTHORS', 'COMMENT', 'CONSRTM', 'DBLINK',
            'DBSOURCE', 'DEFINITION', 'JOURNAL', 'KEYWORDS', 'MEDLINE', 'NID',
            'ORGANISM', 'PID', 'PROJECT', 'PUBMED',
            'REMARK', 'SEGMENT', 'SOURCE', 'TITLE', 'VERSION',
            'VERSION;gi', 'LOCUS;date']
    feature_keys, qual_keys = parse_feature_table(ft_file)

    # -- terms to map them two, along with dictionaries for the existing or
    # hand defined mappings
    dc_terms = parse_dc_terms(dc_file)
    dc_ontology = OntologyGroup("http://purl.org/dc/terms/", dc_terms)
    so_terms = parse_so_terms(so_file)
    so_ontology = OntologyGroup("http://purl.org/obo/owl/SO", so_terms)
    so_ft_map = parse_so_ft_map(so_ft_map_file)
    so_ontology.add_map('feature', 
            'http://www.sequenceontology.org/mappings/FT_SO_map.txt', so_ft_map)
    rdfs_terms = parse_rdfs_terms(rdfs_file)
    rdfs_ontology = OntologyGroup('http://www.w3.org/2000/01/rdf-schema#',
            rdfs_terms)
    me_ft_map = dict(
            rep_origin = 'origin_of_replication',
            unsure = 'sequence_uncertainty',
            conflict = 'sequence_conflict',
            GC_signal = 'GC_rich_promoter_region',
            mat_peptide = 'mature_protein_region',
            C_region = 'C_cluster',
            J_segment = 'J_gene',
            N_region = '',
            S_region = '',
            V_region = 'V_cluster',
            V_segment = 'V_gene'
            )
    me_ft_dc_map = dict(
            old_sequence = 'replaces')
    me_hd_map = {'ACCESSION': 'databank_entry'}
    me_hd_dc_map = {'DEFINITION' : 'description',
                    'VERSION' : 'hasVersion',
                    'VERSION;gi' : 'identifier',
                    'KEYWORDS' : 'subject',
                    'LOCUS;date' : 'created',
                    'PUBMED' : 'relation',
                    'JOURNAL' : 'source',
                    'AUTHORS' : 'contributor',
                    'CONSRTM' : 'creator',
                    'ORGANISM' : '',
                    'SEGMENT' : 'isPartOf',
                    'DBLINK' : 'relation',
                    'DBSOURCE' : 'relation',
                    'MEDLINE' : 'relation',
                    'NID' : 'relation',
                    'PID' : 'relation',
                    'PROJECT' : 'relation'}
    me_ql_map = dict(
            bio_material = 'biomaterial_region',
            bound_moiety = 'bound_by_factor',
            codon_start = 'coding_start',
            direction = 'direction_attribute',
            experiment = 'experimental_result_region',
            macronuclear = 'macronuclear_sequence',
            map = 'fragment_assembly',
            mobile_element = 'integrated_mobile_genetic_element',
            mod_base = 'modified_base_site',
            mol_type = 'sequence_attribute',
            ncRNA_class = 'ncRNA',
            organelle = 'organelle_sequence',
            PCR_primers = 'primer',
            proviral = 'proviral_region',
            pseudo = 'pseudogene',
            rearranged = 'rearranged_at_DNA_level',
            satellite = 'satellite_DNA',
            segment = 'gene_segment',
            rpt_family = 'repeat_family',
            rpt_type = 'repeat_unit',
            rpt_unit_range = 'repeat_region',
            rpt_unit_seq = 'repeat_component',
            tag_peptide = 'cleaved_peptide_region',
            trans_splicing = 'trans_spliced',
            translation = 'polypeptide',
            )
    me_ql_dc_map = dict(
            citation = 'bibliographicCitation',
            gene_synonym = 'alternative',
            identified_by = 'creator',
            label = 'alternative',
            locus_tag = 'alternative',
            old_locus_tag = 'replaces',
            replace = 'isReplacedBy',
            db_xref = 'relation',
            compare = 'relation',
            EC_number = 'relation',
            protein_id = 'relation',
            product = 'alternative',
            standard_name = 'alternative',
            number = 'coverage',
            function = 'description',
            )
    me_ql_rdfs_map = dict(
            note = 'comment',
            )
    ql_make_no_sense_list = ['number']
    so_ontology.add_map('header', 'Brad', me_hd_map)
    so_ontology.add_map('feature', 'Brad', me_ft_map)
    so_ontology.add_map('qualifier', 'Brad', me_ql_map)
    dc_ontology.add_map('header', 'Brad', me_hd_dc_map)
    dc_ontology.add_map('feature', 'Brad', me_ft_dc_map)
    dc_ontology.add_map('qualifier', 'Brad', me_ql_dc_map)
    rdfs_ontology.add_map('qualifier', 'Brad', me_ql_rdfs_map)

    # -- write out the mappings in each of the categories
    with open("genbank_ontology_map.txt", "w") as out_handle:
        out_writer = csv.writer(out_handle, delimiter="\t")
        out_writer.writerow(['gb section', 'identifier', 'ontology',
            'namespace', 'evidence'])
        match_keys_to_ontology('header', header_keys, [so_ontology, dc_ontology,
            rdfs_ontology], out_writer)
        match_keys_to_ontology('feature', feature_keys, [so_ontology,
            dc_ontology, rdfs_ontology], out_writer)
        match_keys_to_ontology('qualifier', qual_keys, [so_ontology,
            dc_ontology, rdfs_ontology], out_writer)

class OntologyGroup:
    def __init__(self, namespace, terms):
        self.ns = namespace
        self.terms = terms
        self._maps = collections.defaultdict(lambda: [])

    def add_map(self, key_type, origin, key_map):
        """Add a mapping of keys to terms within this ontology.
        """
        self._maps[key_type].append((origin, key_map))

    def normalized_terms(self):
        """Retrieve the terms all lower cased and with extra items removed.
        """
        lower_so_terms = {}
        for term in self.terms:
            lower_so_terms[self._normal_term(term)] = term
        return lower_so_terms

    def _normal_term(self, term):
        return term.replace("_", "").lower() 

    def match_key_to_ontology(self, key_type, cur_key):
        normal_terms = self.normalized_terms()
        # try to get it from a dictionary
        for map_origin, check_map in self._maps[key_type]:
            try:
                match_key = check_map[cur_key]
            except KeyError:
                match_key = None
            if (match_key and
                    self._normal_term(match_key) in normal_terms.keys()):
                return match_key, map_origin
        # try to match it by name
        if self._normal_term(cur_key) in normal_terms.keys():
            match_key = normal_terms[self._normal_term(cur_key)]
            return match_key, 'name_match'
        #in_terms = [x for x in normal_terms.keys() if
        #        x.find(self._normal_term(cur_key)) != -1]
        #if len(in_terms) > 0:
        #    print '***', cur_key, in_terms, self.ns
        # could not find a match
        return None, ''

def match_keys_to_ontology(key_type, keys, ontologies, out_writer):
    no_matches = []
    for cur_key in keys:
        found_match = False
        for ontology in ontologies:
            match_key, origin = ontology.match_key_to_ontology(key_type, cur_key)
            if match_key is not None:
                out_writer.writerow([key_type, cur_key, match_key, ontology.ns,
                    origin])
                found_match = True
                break
        if not found_match:
            no_matches.append(cur_key)
    for no_key in no_matches:
        out_writer.writerow([key_type, no_key])

def _parse_terms_from_rdf(in_file, prefix):
    graph = Graph()
    graph.parse(in_file)
    sparta_store = sparta.ThingFactory(graph)
    subjects = []
    for subj in graph.subjects():
        if subj not in subjects:
            subjects.append(subj)
    terms = []
    for subj in subjects:
        rdf_item = sparta_store(subj)
        str_item = str(rdf_item.get_id().replace(prefix + "_", ""))
        if str_item:
            terms.append(str_item)
    terms.sort()
    #print terms
    return terms

def parse_rdfs_terms(rdfs_file):
    return _parse_terms_from_rdf(rdfs_file, "rdfs")

def parse_dc_terms(dc_file):
    """Retrieve a list of Dublin core terms from the RDF file.
    """
    return _parse_terms_from_rdf(dc_file, "dcterms")

def parse_so_terms(so_file):
    """Retrieve all available Sequence Ontology terms from the file.
    """
    so_terms = []
    with open(so_file) as in_handle:
        for line in in_handle:
            if line.find('name:') == 0:
                name = line[5:].strip()
                so_terms.append(name)
    return so_terms

def parse_so_ft_map(so_ft_map_file):
    """Parse out mappings between feature keys and SO.
    """
    so_ft_map = {}
    with open(so_ft_map_file) as in_handle:
        in_handle.readline()
        for line in in_handle:
            parts = line.split()
            if parts[1] not in ['undefined']:
                so_ft_map[parts[0]] = parts[1]
    return so_ft_map

def parse_feature_table(ft_file):
    """Parse all available features and qualifiers from the FT definition.

    This is ugly and parses it straight out of the HTML but this is much easier
    than trying to get it from the specs.
    """
    feature_keys = []
    qual_keys = []
    with open(ft_file) as ft_handle:
        in_feature_region = False
        for line in ft_handle:
            if in_feature_region:
                if line.strip() == "":
                    in_feature_region = False
                else:
                    qual_key, feature_key = line.strip().split()
                    qual_keys.append(qual_key)
                    feature_keys.append(feature_key)
            elif line.find('QUALIFIER FEATURE KEY') == 0:
                in_feature_region = True
    qual_keys = list(set(qual_keys))
    qual_keys = [k.replace('/', '') for k in qual_keys]
    feature_keys = list(set(feature_keys))
    qual_keys.sort()
    feature_keys.sort()
    return feature_keys, qual_keys

if __name__ == "__main__":
    ft_file = os.path.join("feature_table", "FT_index.html")
    so_ft_map_file = os.path.join("feature_table", "FT_SO_map.txt")
    so_file = os.path.join("ontologies", "so.obo")
    dc_file = os.path.join("ontologies", "dcterms.rdf")
    rdfs_file = os.path.join("ontologies", "rdf-schema.rdf")
    main(ft_file, so_ft_map_file, so_file, dc_file, rdfs_file)
