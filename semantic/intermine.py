"""Provide access to servers running Intermine web services.

http://www.intermine.org/wiki/WebService

Intermine is a open source database used for holding experimental data from a
number of model organisms. For instance, the modENCODE project makes their data
available at modMine:

http://intermine.modencode.org/

Queries to do:

    - by lab, affiliation, PI name
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
        #print query
        req = urllib2.Request(self._base,
                urllib.urlencode(dict(query=query)))
        response = urllib2.urlopen(req)
        vals = []
        for line in response:
            parts = line.split('\t')
            parts[-1] = parts[-1].strip()
            vals.append(parts)
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
                et.SubElement(node, "constraint",
                        op=builder.get_compare_op(fname),
                        value=fval, code=name)
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

class _AbstractBuilder:
    """Base class to derive specific query builders from.
    """
    def __init__(self, paths):
        """Provide an initial set of standard items of interest.
        
        paths is a dictionary of object paths to the base class
        of various items. Common retrieval items are automatically added
        for retrieval and selection.
        """
        self._paths = paths
        self._names = {
            "submission_id" : self._path("submission", "DCCid"),
            "organism": self._path("organism", "name"),
            "submission_title" : self._path("submission", "title"),
            "experiment_name" : self._path("experiment", "name"),
            }
        self._back_map = None

        self.attributes = []
        self.filters = []
        
    def _get_back_map(self):
        if self._back_map is None:
            self._back_map = {}
            for key, val in self._names.items():
                self._back_map[val] = key
        return self._back_map

    def get_compare_op(self, out_name):
        """Define comparison operations for different types of values.

        This contains useful operations for various data types.
        """
        back_map = self._get_back_map()
        name = back_map[out_name]
        if name == "start":
            return ">"
        elif name == "end":
            return "<"
        else:
            return "CONTAINS"

    def _path(self, name, attribute):
        return "%s%s" % (self._paths[name], attribute)
    
    def get_out_names(self):
        back_map = self._get_back_map()
        return [back_map[n] for n in self.attributes]

    def available_attributes(self):
        return self._names.keys()

    def add_attributes(self, names):
        if not isinstance(names, (list, tuple)):
            names = [names]
        for name in names:
            self.attributes.append(self._names[name])
    
    def add_filter(self, *args):
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            name_vals = args[0]
        elif len(args) == 2:
            name_vals = [args]
        else:
            raise ValueError("Need a name and value or list of name values")
        filter_group = []
        for name, val in name_vals:
            # A ':' in the name indicates an optional attribute, while a '.'
            # makes it required. All select fields are required to be present
            # since we are selecting on them.
            intermine_name = self._names[name]
            self._names[name] = intermine_name.replace(":", ".")
            if intermine_name.find(":") >= 0:
                #to_swap = ".".join(intermine_name.split(".")[:-1])
                to_swap = intermine_name.split(".")[0]
                new_swap = to_swap.replace(":", ".")
                # change the select in all of our default names
                for sname, sval in self._names.items():
                    if sval.startswith(to_swap):
                        new_val = sval.replace(to_swap, new_swap)
                        self._names[sname] = new_val
                # also swap it in anything we've added
                for i, attr in enumerate(self.attributes):
                    if attr.startswith(to_swap):
                        self.attributes[i] = attr.replace(to_swap, new_swap)
            filter_group.append((self._names[name], val))
        self.filters.append(filter_group)

class LocationQueryBuilder(_AbstractBuilder):
    """Retrieve data associated with a chromosomal region.
    """
    def __init__(self):
        paths = {
                "submission" : "LocatedSequenceFeature:submissions.",
                "organism" : "LocatedSequenceFeature:organism.",
                "experiment": "LocatedSequenceFeature:submissions.experiment",
                }
        _AbstractBuilder.__init__(self, paths)
        self._names.update({
            "chromosome" : "LocatedSequenceFeature:chromosome.name",
            "start" : "LocatedSequenceFeature:chromosomeLocation.start",
            "end": "LocatedSequenceFeature:chromosomeLocation.end",
            "strand": "LocatedSequenceFeature:chromosomeLocation.strand",
        })

class SubmissionQueryBuilder(_AbstractBuilder):
    """Retrieve submissions based on specific submission properties.
    """
    def __init__(self):
        paths = {
                "submission" : "Submission:",
                "organism" : "Submission:organism.",
                "experiment": "Submission:experiment",
                }
        _AbstractBuilder.__init__(self, paths)
        self._names.update({
            "antibody_name" : self._path("submission", "antibodies.name"),
            "cell_line" : self._path("submission", "cellLines.name"),
            "developmental_stage": self._path("submission",
                "developmentalStages.name"),
        })

class ExperimentQueryBuilder(_AbstractBuilder):
    """Provide a mechanism to build high level Experiment queries.
        
    This uses Experiment as a base to perform queries against Itermine.
    The general notion is high level experiment discovery.
    """
    def __init__(self):
        paths = {
                "submission" : "Experiment:project.submissions.",
                "organism" : "Experiment:project.organisms.",
                "experiment" : "Experiment.",
                }
        _AbstractBuilder.__init__(self, paths)
        self._names.update({
          "project_name" : "Experiment:project.name",
          "submission_description": self._path("submission", "description"),
          "experiment_type": self._path("submission", "experimentType"),
        })

    def free_text_filter(self, search_val):
        """Provide a free text style search for the given value.
        """
        self.add_filter([("submission_description", search_val),
                         ("experiment_type", search_val),
                         ("submission_title", search_val)])

# --- Test Code

# Some example XML, for testing
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

    def test_filter_query(self):
        """Provide experiment filtering based on organism and keywords.
        """
        builder = ExperimentQueryBuilder()
        builder.add_attributes([
            "submission_id", "experiment_name"])
        builder.add_filter("organism", "Caenorhabditis elegans")
        builder.free_text_filter("ChIP-seq")

        table = self._server.search(builder)
        print table.dtype.names
        print table
        result = table[0]
        print result['submission_id'], result['experiment_name']

    def test_submission_query(self):
        """Retrieve submissions based on various details of the submission.
        """
        builder = SubmissionQueryBuilder()
        builder.add_attributes(["submission_id", 
            "submission_title", "developmental_stage"])
        builder.add_filter("organism", "Caenorhabditis elegans")
        builder.add_filter("antibody_name", "H3K4me3")
        
        table = self._server.search(builder)
        print table.dtype.names
        print table

    def test_location_query(self):
        """Retrieve submissions with data in particular chromosome locations.
        """
        builder = LocationQueryBuilder()
        builder.add_attributes(["submission_id",
            "submission_title"])
        builder.add_filter("organism", "Caenorhabditis elegans")
        builder.add_filter("chromosome", "I")
        builder.add_filter("start", "5000")
        builder.add_filter("end", "20000")
        
        table = self._server.search(builder)
        print table.dtype.names
        print table
