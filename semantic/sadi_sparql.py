import json
import time
import urllib, urllib2
from rdflib.Graph import Graph

result_url_base = "http://biordf.net/tmp/"
query_url = "http://biordf.net/cardioSHARE/query"

query = """
PREFIX pred: <http://sadiframework.org/ontologies/predicates.owl#> 
PREFIX uniprot: <http://lsrn.org/UniProt:> 
SELECT ?name WHERE { 
    uniprot:P15923 pred:hasName ?name 
    }
"""

query = " ".join(query.split("\n"))
req = urllib2.Request(query_url, urllib.urlencode(dict(query=query)))
response = urllib2.urlopen(req)

info = json.loads(response.read())
poll_url = query_url + "?" + urllib.urlencode(dict(poll=info["taskId"]))
while 1:
    response = urllib2.urlopen(poll_url)
    poll_text = response.read()
    # got our JSON response -- means we are ready to retrieve
    if poll_text.startswith("{"):
        break
    time.sleep(3)
poll_info = json.loads(poll_text)

results_url = result_url_base + info["taskId"]

g = Graph()
g.parse(results_url)
for stmt in g:
    print stmt
