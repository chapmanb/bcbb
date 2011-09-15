"""Python API to query BioMart SPARQL endpoints.

BioMart SPARQL documentation is at: http://www.biomart.org/rc6_documentation.pdf (p83-85)

Some useful BioMart servers:
BMC:  http://bm-test.res.oicr.on.ca:9084/
ICGC: http://bm-test.res.oicr.on.ca:9085/
"""
import unittest
from xml.etree import ElementTree as ET
from collections import namedtuple

import SPARQLWrapper

class SematicBioMart:
    """Given SPARQL query, retrieve results from remote BioMart SPARQL endpoint.
    """
    def __init__(self, base_url):
        self._base_url = base_url
        self._query_url = "http://{url}/martsemantics/{ap}/SPARQLXML/get/"
        self._result_ns = "{http://www.w3.org/2005/sparql-results#}"
        self._prefixes = [
          ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"),
          ("rdfs", "http://www.w3.org/2000/01/rdf-schema#"),
          ("owl", "http://www.w3.org/2002/07/owl#"),
          ("accesspoint", "http://{url}/martsemantics/{ap}/ontology#"),
          ("class", "biomart://{url}/martsemantics/{ap}/ontology/class#"),
          ("dataset", "biomart://{url}/martsemantics/{ap}/ontology/dataset#"),
          ("attribute", "biomart://{url}/martsemantics/{ap}/ontology/attribute#"),
         ]

    def _get_prefix_str(self, access_point):
        """Convert to prefix strings: PREFIX owl: <http://www.w3.org/2002/07/owl#>
        """
        out = ["PREFIX {name}: <{url}>".format(name=n, url=u.format(url=self._base_url,
                                                                    ap=access_point))
               for n, u in self._prefixes]
        return "\n".join(out) + "\n\n"

    def _do_query(self, statement, access_point):
        sparql = SPARQLWrapper.SPARQLWrapper(self._query_url.format(url=self._base_url,
                                                                    ap=access_point))
        sparql.setReturnFormat(SPARQLWrapper.XML)
        query = self._get_prefix_str(access_point) + statement
        #print query
        sparql.setQuery(query)
        #for line in sparql.query().response:
        #    print line.rstrip()
        return self._parse_response(sparql.query().response)

    def _parse_response(self, in_handle):
        tree = ET.parse(in_handle).getroot()
        out = []
        Result = None
        for results in tree.findall(self._result_ns + "results"):
            for result in results.findall(self._result_ns + "result"):
                names = []
                vals = []
                for bind in result.findall(self._result_ns + "binding"):
                    names.append(self._parse_biomart_name(bind.get("name")))
                    vals.append(bind.find(self._result_ns + "literal").text)
                if Result is None:
                    Result = namedtuple("Result", names)
                out.append(Result(*vals))
        return out 

    def _parse_biomart_name(self, name):
        """Retrieve last part of long BioMart name as attribute.
        """
        return name.split("__")[-1]

    def search(self, builder):
        access_point = "{dataset}_config".format(dataset=builder.dataset)
        return self._do_query(builder.sparql(), access_point)

class BioMartQueryBuilder:
    def __init__(self, dataset, sub_dataset):
        self.dataset = dataset
        self.sub_dataset = sub_dataset
        self._attrs = []
        self._filters = []
        mart_content = {
            ("simple_somatic_mutation", "dm") : ["consequence_type", "validation_status",
                                                 "aa_mutation", "gene_affected",
                                                 "probability", "mutation"],
            ("feature", "main") : ["chromosome", "chromosome_start", "chromosome_end"]}
        self._attr_lookup = self._content_lookup(mart_content)

    def _content_lookup(self, mart_content):
        """Provide mapping of attributes to content and type values.
        """
        out = {}
        for (content, mart_type), attrs in mart_content.iteritems():
            for attr in attrs:
                out[attr] = (content, mart_type)
        return out

    def add_attributes(self, attrs):
        self._attrs.extend([self._biomart_name(a) for a in attrs])

    def add_filter(self, attr, val):
        self._filters.append((self._biomart_name(attr), val))

    def _biomart_name(self, attr):
        """Generate long BioMart attribute name from shortened name.

        These consist of 4 parts (section 4.1 http://www.biomart.org/install.html):
        - dataset
        - content: a free text description
        - type: main or dm (dimension tables)
        - short attribute name
        """
        content, mart_type = self._attr_lookup[attr]
        return "__".join([self.dataset, content, mart_type, attr])

    def sparql(self):
        """Retrieve the SPARQL query for currently set attributes and filters.
        """
        return "\n".join([self._sparql_select(),
                          self._sparql_from(),
                          self._sparql_where()])

    def _sparql_select(self):
        """Build SELECT portion of SPARQL query.
        """
        return "SELECT {attrs}".format(attrs=" ".join(["?{0}".format(a) for a in self._attrs]))

    def _sparql_from(self):
        """Build FROM portion of SPARQL query
        """
        return "FROM dataset:{main}_{sub}".format(main=self.dataset, sub=self.sub_dataset)

    def _sparql_where(self):
        """WHERE clause of SPARQL query.
        """
        select_lines = []
        for attr, val in self._filters:
            select_lines.append('?x attribute:{attr} "{value}" .'.format(attr=attr, value=val))
        for attr in self._attrs:
            select_lines.append("?x attribute:{attr} ?{attr} .".format(attr=attr))
        return "WHERE {{\n{0}\n}}".format("\n".join(select_lines))

class BioMartSematicTest(unittest.TestCase):
    def test_snp_positions(self):
        """Find positions and changes of validated SNPs
        """
        builder = BioMartQueryBuilder("snp", "jpNCCLiver")
        builder.add_attributes(["chromosome", "chromosome_start", "chromosome_end",
                                "aa_mutation", "gene_affected", "probability", "mutation"])
        builder.add_filter("consequence_type", "non_synonymous_coding")
        builder.add_filter("validation_status", "validated")
        icgc_server = SematicBioMart("bm-test.res.oicr.on.ca:9085")
        results = icgc_server.search(builder)
        print results[0]
        assert results[0].chromosome == "1"
        assert results[0].aa_mutation == "D>Y"

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
