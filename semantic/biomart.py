"""Python API to query BioMart SPARQL endpoints.

BioMart SPARQL documentation is at: http://www.biomart.org/rc6_documentation.pdf (p83-85)

Some useful BioMart servers:
BMC:  http://bm-test.res.oicr.on.ca:9084/
ICGC: http://bm-test.res.oicr.on.ca:9085/
"""
import unittest

from SPARQLWrapper import SPARQLWrapper

class SematicBioMart:
    """Given SPARQL query, retrieve results from remote BioMart SPARQL endpoint.
    """
    def __init__(self, base_url):
        self._base_url = base_url
        self._access_point = "snp_config"
        self._query_url = "{url}/martsemantics/{ap}/SPARQLXML/get/".format(
            url=base_url, ap=self._access_point)

        self._prefixes = [
          ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"),
          ("rdfs", "http://www.w3.org/2000/01/rdf-schema#"),
          ("owl","http://www.w3.org/2002/07/owl#"),
          ("accesspoint", "http://{url}/martsemantics/{ap}/ontology#"),
          ("class","biomart://{url}/martsemantics/{ap}/ontology/class#"),
          ("dataset","biomart://{url}/martsemantics/{ap}/ontology/dataset#"),
          ("attribute","biomart://{url}/martsemantics/{ap}/ontology/attribute#"),
          ]

    def _get_prefix_str(self):
        """Convert to prefix strings: PREFIX owl: <http://www.w3.org/2002/07/owl#>
        """
        out = ["PREFIX {name}: <{url}>".format(name=n, url=u.format(url=self._base_url,
                                                                    ap=self._access_point))
               for n, u in self._prefixes]
        return "\n".join(out) + "\n\n"

    def do_query(self, statement):
        sparql = SPARQLWrapper(self._query_url)
        query = self._get_prefix_str() + statement
        print query
        sparql.setQuery(query)
        print sparql.query().geturl()
        for info in sparql.query():
            print info

class BioMartQueryBuilder:
    def __init__(self, dataset):
        pass

    def add_attributes(self, attrs):
        pass

    def add_filter(self, attr, val):
        pass

example_query = """
SELECT ?snp__feature__main__chromosome ?snp__feature__main__chromosome_start ?snp__feature__main__chromosome_end
       ?snp__simple_somatic_mutation__dm__aa_mutation ?snp__simple_somatic_mutation__dm__gene_affected
       ?snp__simple_somatic_mutation__dm__probability ?snp__simple_somatic_mutation__dm__mutation
FROM dataset:snp_jpNCCLiver
WHERE {
  ?x attribute:snp__simple_somatic_mutation__dm__consequence_type "non_synonymous_coding" .
  ?x attribute:snp__simple_somatic_mutation__dm__validation_status "validated" .
  ?x attribute:snp__feature__main__chromosome ?snp__feature__main__chromosome .
  ?x attribute:snp__feature__main__chromosome_start ?snp__feature__main__chromosome_start .
  ?x attribute:snp__feature__main__chromosome_end ?snp__feature__main__chromosome_end .
  ?x attribute:snp__simple_somatic_mutation__dm__aa_mutation ?snp__simple_somatic_mutation__dm__aa_mutation .
  ?x attribute:snp__simple_somatic_mutation__dm__gene_affected ?snp__simple_somatic_mutation__dm__gene_affected .
  ?x attribute:snp__simple_somatic_mutation__dm__probability ?snp__simple_somatic_mutation__dm__probability .
  ?x attribute:snp__simple_somatic_mutation__dm__mutation ?snp__simple_somatic_mutation__dm__mutation
}
"""

icgc_server = SematicBioMart("http://bm-test.res.oicr.on.ca:9085")
icgc_server.do_query(example_query)
