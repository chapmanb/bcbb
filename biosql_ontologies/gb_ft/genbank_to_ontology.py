"""Provide a (semi) automated mapping between GenBank and ontologies.

The goal is to provide an ontology namespace and term for each item
of a standard GenBank file.
"""
from __future__ import with_statement

def main(ft_file, so_ft_map_file, so_file):
    # pulled by hand from the Biopython parser
    header_keys = ['ACCESSION', 'AUTHORS', 'COMMENT', 'CONSRTM', 'DBLINK',
            'DBSOURCE', 'DEFINITION', 'JOURNAL', 'KEYWORDS', 'MEDLINE', 'NID',
            'ORGANISM', 'ORIGIN', 'PID', 'PROJECT', 'PUBMED', 'REFERENCE',
            'REMARK', 'SEGMENT', 'SOURCE', 'TITLE', 'VERSION']
    so_terms = parse_so_terms(so_file)
    feature_keys, qual_keys = parse_feature_table(ft_file)
    me_ft_map = dict(
            rep_origin = 'origin_of_replication',
            unsure = 'sequence_uncertainty',
            conflict = 'sequence_conflict',
            GC_signal = 'GC_rich_promoter_region',
            mat_peptide = 'mature_protein_region',
            )
    so_ft_map = parse_so_ft_map(so_ft_map_file)
    so_ft_map.update(me_ft_map)
    match_keys_to_ontology('header', header_keys, so_terms, {})
    match_keys_to_ontology('feature', feature_keys, so_terms, so_ft_map)
    match_keys_to_ontology('qualifier', qual_keys, so_terms, {})

def match_keys_to_ontology(key_type, keys, so_terms, pre_map):
    lower_so_terms = [t.lower() for t in so_terms]
    print '***', key_type
    nos = []
    for key in keys:
        try:
            pre_key = pre_map[key]
            if pre_key in so_terms:
                print 'pre', key
            else:
                raise KeyError
        except KeyError:
            if key.lower() in lower_so_terms:
                print 'SO ', key
            else:
                nos.append(key)
    print 'no mapping', nos

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
    ft_file = "FT_index.html"
    so_ft_map_file = "FT_SO_map.txt"
    so_file = "so.obo"
    main(ft_file, so_ft_map_file, so_file)
