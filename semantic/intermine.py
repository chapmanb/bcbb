"""Provide access to servers running Intermine web services.

http://www.intermine.org/wiki/WebService

Intermine is a open source database used for holding experimental data from a
number of model organisms. For instance, the modENCODE project makes their data
available at modMine:

http://intermine.modencode.org/
"""
import string
import unittest
import StringIO
import urllib, urllib2
from xml.etree import ElementTree as et

import numpy

class Intermine:
    """Provide query based access to data through Intermine web services.
    """
    def __init__(self, base_url):
        self._base = "%s/query/service/query/results" % base_url

    def _do_query(self, query):
        req = urllib2.Request(self._base,
                urllib.urlencode(dict(query=query)))
        response = urllib2.urlopen(req)
        vals = []
        for line in response:
            vals.append(line.strip().split('\t'))
        return vals

    def search(self, builder):
        """Query intermine and return results based on the provided builder.
        """
        # build our filter statements
        nodes = []
        constraints = []
        i = 0
        for filter_group in builder.filters:
            group_names = []
            for fname, fval in filter_group:
                name = string.uppercase[i]
                group_names.append(name)
                node = et.Element("node", path=fname, type="String")
                et.SubElement(node, "constraint", op="CONTAINS", value=fval,
                              code=name)
                nodes.append(node)
                i += 1
            constraints.append("(%s)" % " or ".join(group_names))
        # now build the query
        query = et.Element('query', model="genomic", 
                view=" ".join(builder.attributes),
                constraintLogic=" and ".join(constraints))
        for node in nodes:
            query.append(node)
        # serialize and send
        query_handle = StringIO.StringIO()
        et.ElementTree(query).write(query_handle)
        vals = self._do_query(query_handle.getvalue())
        term_names = builder.get_out_names()
        return numpy.core.records.array(vals, names=",".join(term_names))

class ExperimentQueryBuilder:
    """Provide a mechanism to build high level Experiment queries.
        
    This uses Experiment as a base to perform queries against Itermine.
    The general notion is high level experiment discovery.
    """
    def __init__(self):
        self._names = {
          "experiment_name" : "Experiment.name",
          "project_name" : "Experiment.project.name",
          "submission_id" : "Experiment.project.submissions.DCCid",
          "submission_description": "Experiment.project.submissions.description",
          "experiment_type": "Experiment.project.submissions.experimentType",
          "submission_title": "Experiment.project.submissions.title",
          "organism_name": "Experiment.project.organisms.name",
        }

        self.attributes = []
        self.filters = []

    def get_out_names(self):
        back_map = {}
        for key, val in self._names.items():
            back_map[val] = key
        return [back_map[n] for n in self.attributes]

    def available_names(self):
        return self._names.keys()

    def add_attributes(self, names):
        if not isinstance(names, (list, tuple)):
            names = [names]
        for name in names:
            self.attributes.append(self._names[name])
    
    def filter_organism(self, org):
        """Filter our query to include the given organism.
        """
        self.add_filter([("organism_name", org)])

    def free_text_filter(self, search_val):
        """Provide a free text style search for the given value.
        """
        self.add_filter([("submission_description", search_val),
                         ("experiment_type", search_val),
                         ("submission_title", search_val)])

    def add_filter(self, name_vals):
        filter_group = []
        for name, val in name_vals:
            filter_group.append((self._names[name], val))
        self.filters.append(filter_group)

q = """
<query name="" model="genomic" view=" Experiment.project.submissions.DCCid Experiment.name Experiment.project.name" constraintLogic="B and (A or C or D)">
  <node path="Experiment.project.organisms.name" type="String">
    <constraint op="=" value="Caenorhabditis elegans" description="" identifier="" code="B" extraValue="">
    </constraint>
  </node>
  <node path="Experiment.project.submissions.description" type="String">
    <constraint op="CONTAINS" value="ChIP-seq" description="" identifier="" code="A" extraValue="">
    </constraint>
  </node>
  <node path="Experiment.project.submissions.experimentType" type="String">
    <constraint op="=" value="ChIP-seq" description="" identifier="" code="C" extraValue="">
    </constraint>
  </node>
  <node path="Experiment.project.submissions.title" type="String">
    <constraint op="CONTAINS" value="ChIP-seq" description="" identifier="" code="D" extraValue="">
    </constraint>
  </node>
</query>
"""

class IntermineTest(unittest.TestCase):
    def setUp(self):
        self._server = Intermine("http://intermine.modencode.org")
       
    def test_query(self):
        """Simple string based query with Intermine XML.
        """
        vals = self._server._do_query(q)
        print vals

    def test_filter_query(self):
        """Provide filtering based on organism and keywords.
        """
        builder = ExperimentQueryBuilder()
        builder.add_attributes([
            "submission_id", "experiment_name", "project_name"])
        builder.filter_organism("Caenorhabditis elegans")
        builder.free_text_filter("ChIP-seq")
        table = self._server.search(builder)
        print table.dtype.names
        result = table[0]
        print result['submission_id'], result['experiment_name']
