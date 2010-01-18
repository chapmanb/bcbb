import urllib
from rdflib import ConjunctiveGraph as Graph
import sparta

url = 'http://www.gopubmed.org/GoMeshPubMed/gomeshpubmed/Search/RDF?q=18463287&type=RdfExportAll'
gopubmed_handle = urllib.urlopen(url)
graph = Graph()
graph.parse(gopubmed_handle)
gopubmed_handle.close()

graph_subjects = list(set(graph.subjects()))
sparta_factory = sparta.ThingFactory(graph)
for subject in graph_subjects:
    sparta_graph = sparta_factory(subject)
    print subject, [unicode(i) for i in sparta_graph.dc_title][0]
