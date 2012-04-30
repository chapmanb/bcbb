"""Provide a client API to do queries on resources at Semantic Systems Biology.

This wraps the SPARQL query endpoint for Biogateway:

http://www.semantic-systems-biology.org/biogateway
"""
import unittest

import numpy
from SPARQLWrapper import SPARQLWrapper, JSON

class Biogateway:
    """Provide a query builder for getting Biogateway resources.
    """
    def __init__(self):
        self._base_url = "http://www.semantic-systems-biology.org/"
        self._query_url = "%s/biogateway/endpoint" % self._base_url

        self._org_map = None

        self._ns = {
                "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
                "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
                "ssb": "%sSSB#" % self._base_url,
                }

    def _query_header(self):
        header = "BASE   <%s>\n" % self._base_url
        for short, nurl in self._ns.items():
            header += "PREFIX %s:<%s>\n" % (short, nurl)
        return header

    def _strip_ns(self, data):
        if data.startswith("http"):
            for rem in self._ns.values():
                data = data.replace(rem, "")
        return data

    def _do_query(self, query):
        """Perform the actual work of doing a SPARQL query and parsing results.
        """
        sparql = SPARQLWrapper(self._query_url)
        query = self._query_header() + query
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        out = []
        for result in results["results"]["bindings"]:
            cur_out = dict()
            for var in results["head"]["vars"]:
                try:
                    data = result[var]["value"]
                    data = self._strip_ns(data)
                except KeyError:
                    data = ""
                cur_out[var] = data
            out.append(cur_out)
        return out

    def get_organisms(self):
        """Retrieve organisms available in the database.
        """
        query = """
        SELECT distinct ?taxon ?graph
        WHERE {
          GRAPH <metaonto> {
            ?graph ssb:about_taxon ?taxon_id.
          }
          GRAPH <ncbi> {
            ?taxon_id rdfs:label   ?taxon.
          }
          FILTER(?graph != <SSB> && ?graph != <GOA> && ?graph != <SSB_tc>).
        }
        """
        if self._org_map is None:
            self._org_map = dict()
            for taxon_info in self._do_query(query):
                if not taxon_info["graph"].endswith("_tc"):
                    if not self._org_map.has_key(taxon_info["taxon"]):
                        self._org_map[taxon_info["taxon"]] = taxon_info["graph"]
        orgs = self._org_map.keys()
        orgs.sort()
        return orgs

    def search(self, builder):
        """Retrieve protein IDs from Uniprot with given selection terms.

        Returns a numpy named record array representing the table of results:

        http://www.scipy.org/RecordArrays
        """
        self.get_organisms() # load up our organism mapping
        #XXX for testing
        #self._org_map = {"Homo sapiens" : "25.H_sapiens"}
        all_terms = builder.attributes + builder.filters
        stmt = self._get_sparql_piece("SELECT distinct", "select",
                all_terms, " ")
        stmt += "\nWHERE {"
        if builder.organism is not None:
            stmt += self._get_sparql_piece("GRAPH <%s> {" %
                    self._org_map[builder.organism],
                    "org_graph", all_terms, "\n", "}\n")
        for (graph_name, attr_name) in [("uniprot_sprot", "uniprot_graph"),
                                        ("SSB", "ssb"),
                                        ("gene_ontology_edit", "go_graph"),
                                        ("evidence_code", "evidence_graph")]:
            stmt += self._get_sparql_piece("GRAPH <%s> {" % graph_name,
                    attr_name, all_terms, "\n", "}\n")
        stmt += self._get_sparql_piece("", "to_filter", all_terms, "\n")
        stmt += "\n}"
        results = self._do_query(stmt)
        term_names = [t.select[1:] for t in all_terms if t.select]
        vals = []
        for r in results:
            vals.append([r[n] for n in term_names])
        if len(vals) > 0:
            vals = numpy.core.records.array(vals, names=",".join(term_names))
        else:
            vals = None
        return vals

    def _get_sparql_piece(self, base, attr, terms, join, end = ""):
        stmt = ""
        for term in terms:
            to_add = getattr(term, attr, None)
            if to_add:
                stmt += "%s%s" % (join, to_add)
        if stmt:
            stmt = base + stmt + end
        return stmt

class _AbstractBuilder:
    """Base class to derive specific query builders from.
    """
    def __init__(self):
        self._terms = {}
        self._selects = {}
        self.organism = None
        self.attributes = []
        self.filters = []

    def available_attributes(self):
        return sorted(self._terms.keys())

    def add_attributes(self, names):
        for n in names:
            self.attributes.append(self._terms[n]())

    def available_filters(self):
        return sorted(self._selects.keys())

    def add_filter(self, *args):
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            name_vals = args[0]
        elif len(args) == 2:
            name_vals = [args]
        else:
            raise ValueError("Need a name and value or list of name values")
        for name, val in name_vals:
            self.filters.append(self._selects[name](val))

class UniProtGOQueryBuilder(_AbstractBuilder):
    """Build queries for retrieval against UniProt and GO.

    Biogateway contains the SwissProt database integrated with associated Gene
    Ontology terms. This provides a builder to help attain common linked
    information of interest.
    """
    def __init__(self, organism):
        _AbstractBuilder.__init__(self)
        for tclass in [_RetrieveProtein, _RetrieveInteractor,
                _RetrieveGeneName]:
            self._terms[tclass.select] = tclass
        for tclass in [_SelectByGOTerm, _SelectByDisease]:
            self._selects[tclass.select] = tclass

        self.organism = organism

class ReferenceBuilder(_AbstractBuilder):
    """Build queries to retrieve GO annotations and references.
    """
    def __init__(self):
        _AbstractBuilder.__init__(self)
        for tclass in [_RetrieveReference, _RetrieveEvidence,
                _RetrieveGODescription]:
            self._terms[tclass.select] = tclass
        for tclass in [_SelectByProteinName]:
            self._selects[tclass.select] = tclass

# -- Useful definitions for retrieving and selecting by common items

FILTER_BASE = "FILTER regex(str(%s), '%s')."

class _RetrieveProtein:
    select = "protein_name"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.org_graph = "?protein_id rdfs:label %s." % self.select

class _RetrieveInteractor:
    select = "interactor"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.org_graph = "?interactor_id rdfs:label %s." % self.select
        self.uniprot_graph = """OPTIONAL {
         ?protein_id ssb:interacts_with ?interactor_id.
        }
        """

class _RetrieveGeneName:
    select = "gene_name"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.org_graph = """
        OPTIONAL{
          ?protein_id ssb:encoded_by %s.
          }
         """ % self.select

class _RetrieveReference:
    select = "reference"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.ssb = """
        ?GOA_triple rdf:subject ?protein_id.
        ?GOA_triple ssb:supported_by ?support_node.
        ?support_node ssb:refer %s.
        """ % self.select

class _RetrieveEvidence:
    select = "evidence"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.ssb = """
        ?GOA_triple rdf:subject ?protein_id.
        ?GOA_triple ssb:supported_by ?support_node.
        ?support_node ssb:has_evidence ?evidence_id.
        """
        self.evidence_graph = """
        ?evidence_id rdfs:label %s.
        ?evidence_id a ?type1.
        """ % self.select

class _RetrieveGODescription:
    select = "go_desc"
    def __init__(self):
        self.select = "?%s" % self.__class__.select
        self.ssb = "?GOA_triple rdf:object ?object."
        self.go_graph = """
          ?object rdfs:label %s.
          ?object a ?type2.
          """ % self.select

class _SelectByGOTerm:
    """Provide selection of information based on GO keywords.
    """
    select = "GO_term"
    def __init__(self, keyword):
        self.select = "?%s" % self.__class__.select
        self.to_filter = FILTER_BASE % (self.select, keyword)
        self.go_graph = "?GO_id rdfs:label %s.\n" % self.select
        self.org_graph = "?protein_id ?relation_id ?GO_id."

class _SelectByDisease:
    select = "disease_description"
    def __init__(self, keyword):
        self.select = "?%s" % self.__class__.select
        self.to_filter = FILTER_BASE % (self.select, keyword)
        self.uniprot_graph = "?protein_id ssb:disease %s." % self.select

class _SelectByProteinName:
    select = "protein_id"
    def __init__(self, keyword):
        self.select = "?%s" % self.__class__.select
        self.to_filter = "FILTER regex(?found_in_name,'%s','i')" % keyword
        self.uniprot_graph = """
          {%s ssb:name ?found_in_name.}
                UNION
          {%s ssb:mnemonic ?found_in_name.}
                UNION
          {%s ssb:encoded_by ?found_in_name.}
        """ % ((self.select,) * 3)

class BiogatewayTest(unittest.TestCase):
    """Test retrieval from the Biogateway Systems Bio server.
    """
    def setUp(self):
        self._server = Biogateway()

    def test_organism(self):
        """Retrieve organisms available for querying.
        """
        orgs = self._server.get_organisms()
        print orgs[:5]

    def test_search_query(self):
        """Build a query for searching based on GO terms and diseases.
        """
        builder = UniProtGOQueryBuilder("Homo sapiens")
        builder.add_attributes(["protein_name", "interactor", "gene_name"])
        builder.add_filter("GO_term", "insulin")
        builder.add_filter("disease_description", "diabetes")
        results = self._server.search(builder)
        print len(results), results.dtype.names
        result = results[0]
        print result['protein_name'], result['GO_term'], result['interactor'], \
              result['disease_description']

    def needsupdatetest_reference_query(self):
        """Retrieve GO associations from a UniProt name.
        """
        builder = ReferenceBuilder()
        builder.add_attributes(["reference"])
        builder.add_filter("protein_id", "1433B_HUMAN")
        results = self._server.search(builder)
        print len(results), results.dtype.names
        result = results[0]
        print result
        print result['protein_id'], result['reference']
