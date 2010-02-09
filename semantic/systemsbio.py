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

    def _query_header(self):
        return """
          BASE   <%s>
          PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
          PREFIX ssb:<%sSSB#>
        """ % ((self._base_url,) * 2)

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
                    self._org_map[taxon_info["taxon"]] = taxon_info["graph"]
        orgs = self._org_map.keys()
        orgs.sort()
        return orgs

    def search(self, organism, to_get, select_terms):
        """Retrieve protein IDs from Uniprot with given selection terms.

        Returns a numpy named record array representing the table of results:

        http://www.scipy.org/RecordArrays
        """
        self.get_organisms() # load up our organism mapping
        #XXX for testing
        #self._org_map = {"Homo sapiens" : "25.H_sapiens"}
        all_terms = to_get + select_terms
        stmt = self._get_sparql_piece("SELECT distinct", "select",
                all_terms, " ")
        stmt += "\nWHERE {"
        stmt += self._get_sparql_piece("GRAPH <gene_ontology_edit> {",
                "go_graph", all_terms, "\n", "}\n")
        stmt += self._get_sparql_piece("GRAPH <%s> {" % self._org_map[organism],
                "org_graph", all_terms, "\n", "}\n")
        stmt += self._get_sparql_piece("GRAPH <uniprot_sprot> {",
                "uniprot_graph", all_terms, "\n", "}\n")
        stmt += self._get_sparql_piece("", "to_filter", all_terms, "\n")
        stmt += "\n}"
        results = self._do_query(stmt)
        term_names = [t.select[1:] for t in all_terms if t.select]
        vals = []
        for r in results:
            vals.append([r[n] for n in term_names])
        vals = numpy.core.records.array(vals, names=",".join(term_names))
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

FILTER_BASE = "FILTER regex(str(%s), '%s')."

# -- Useful definitions for retrieving and selecting by common items

class RetrieveProtein:
    def __init__(self):
        self.select = "?protein"
        self.org_graph = "?protein_id rdfs:label %s." % self.select

class RetrieveInteractor:
    def __init__(self):
        self.select = "?interactor"
        self.org_graph = "?interactor_id rdfs:label %s." % self.select
        self.uniprot_graph = """OPTIONAL {
         ?protein_id ssb:interacts_with ?interactor_id.
        }
        """

class SelectByGOTerm:
    """Provide selection of information based on GO keywords.
    """
    def __init__(self, keyword):
        self.select = "?GO_term"
        self.to_filter = FILTER_BASE % (self.select, keyword)
        self.go_graph = "?GO_id rdfs:label %s.\n" % self.select
        self.org_graph = "?protein_id ?relation_id ?GO_id."

class SelectByDisease:
    def __init__(self, keyword):
        self.select = "?disease_description"
        self.to_filter = FILTER_BASE % (self.select, keyword)
        self.uniprot_graph = "?protein_id ssb:disease %s." % self.select

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

    def test_select_query(self):
        """Generic select and filter queries on a database.
        """
        results = self._server.search("Homo sapiens", 
                [RetrieveProtein(), RetrieveInteractor()],
                [SelectByGOTerm('insulin'), SelectByDisease('diabetes')])
        print results.dtype.names
        result = results[0]
        print result['protein'], result['interactor'], result['GO_term']
