#!/usr/bin/env python
"""Prepare a tree representation of
"""
from __future__ import with_statement
import os
import sys
import re
import glob
import shelve
import networkx
import collections
import pygraphviz as pgv

from BeautifulSoup import BeautifulSoup

branch_includes = ["Metazoa", "Nematoda", "Vertebrata", "Mammalia",
        "Bovidae", "Rodentia", "Primates", "Teleostei", "Amphibia",
        "Aves", "Arthropoda", "Cephalochordata", "Cnidaria", "Chordata"]

org_includes = [
        'Homo sapiens',
        'Mus musculus',
        'Rattus norvegicus',
        'Bos taurus',
        'Gallus gallus',
        'Xenopus laevis',
        'Danio rerio',
        'Tetraodon nigroviridis',
        'Ciona savignyi',
        'Ciona intestinalis',
        'Branchiostoma floridae',
        'Caenorhabditis elegans',
        'Anopheles gambiae',
        'Drosophila melanogaster',
        'Nematostella vectensis'
]

# hard code items that have tandem duplications, so we don't have to
# rebuild this here
with_tandem = ["Mus musculus", "Homo sapiens", "Drosophila melanogaster",
        "Danio rerio"]#, "Bos taurus"]

def main(base_dir):
    db_dir = os.path.join(os.getcwd(), "db")
    cur_dbs = get_available_dbs(db_dir)
    #base_dir = "/usr/home/chapmanb/mgh/dan_grau/analysis_domains/"\
    #        "no_interactions_cluster/"
    #base_dir = "/usr/home/chapmanb/mgh/src/wormbase_lite/web/WormBaseLite/"\
    #        "wormbaselite/public/dan/"
    org_classifications = get_organism_classifications(cur_dbs, base_dir)

    gv_graph = build_tax_graph(cur_dbs, org_includes, org_classifications)
    #gv_graph = networkx.to_agraph(tax_graph)
    gv_graph = set_org_layout(gv_graph, org_includes)
    gv_graph.graph_attr["rankdir"] = "LR"
    gv_graph.graph_attr["nodesep"] = ".1"
    gv_graph.graph_attr["ranksep"] = ".1"
    gv_graph.graph_attr["ordering"] = "in"
    gv_graph.node_attr["fontname"] = "courier"
    #gv_graph.node_attr["fontname"] = "Helvetica"
    gv_graph.edge_attr["arrowhead"] = "none"
    gv_graph.write("org_charge.dot")
    gv_graph.layout(prog='dot')
    #gv_graph.layout(prog='neato')
    gv_graph.draw(os.path.join(base_dir, 'tax_results_summary.png'))

def get_organism_classifications(cur_dbs, base_dir):
    """Determine organism classifications based on grouping of genes in clusters

    This groups proteins from an organism into a 4x4 matrix to summarize the
    likely activity for that organism. The groupings are based on k-means
    clustering of genes into the matrix:

        Chromo         |  Chromo
        Cbx2 (active)  |  Pc/Cbx7 (not-active)
        -------------------------------------
        Zinc ring      |  Zing ring
        Psc (active)   |  Bmi1 (not-active)
    """
    file_base = os.path.join(base_dir, "*-cluster-%s.html")
    groups = [(("chromo", "active"), ["Cbx2"]),
              (("chromo", "non-active"), ["Pc"]),#, "Cbx7"]),
              (("zinc-ring", "active"), ["Psc"]),#, "Suz2"]),
              (("zinc-ring", "non-active"), ["Bmi1"])]
    org_info = collections.defaultdict(lambda: collections.defaultdict(
        lambda: []))
    for cur_group, group_names in groups:
        all_ids = []
        for group_name in group_names:
            all_ids.extend(ids_from_group(group_name, file_base))
        print cur_group, len(all_ids)
        for cur_id in all_ids:
            db_rec = get_db_rec(cur_id, cur_dbs)
            org_info[db_rec["org_scientific_name"]][cur_group].append(cur_id)
    org_classifiers = {}
    for org in org_info.keys():
        total_count = sum([len(v) for (d, a), v in org_info[org].items()])
        active_count = sum([len(v) for (d, a), v in org_info[org].items() 
            if a == "active"])
        if active_count > 0 or (total_count - active_count) > 1:
            print org, total_count
            group_info = {}
            for cur_group, group_names in groups:
                items = org_info[org][cur_group]
                group_info[cur_group] = items
                print '\t', cur_group, len(items)
            org_classifiers[org] = classify_org(org, group_info)
    return org_classifiers

def classify_org(org, group_info):
    """Given counts of items that fall into various groups of our matrix.

    This provides the final classification work that can go into a display
    of the results, and summarizes three details:

    - Chromo-domain proteins (active or not active or not present)
    - Zinc-ring proteins (active or not active or not present)
    - Tandem duplications (identified or not)
    """
    ch_labels = {("chromo", "active") : 'ch',
                 ("chromo", "non-active") : '-',
                 ("zinc-ring", "active"): 'zn',
                 ("zinc-ring", "non-active"): '-'}
    all_groups = [name for name, items in group_info.items() if len(items) > 0]
    classifiers = []
    for base_name, short_name in (("chromo", "ch"), ("zinc-ring", "zn")):
        if (base_name, "active") in all_groups:
            new_class = '%s (%s)' % (short_name,
                len(group_info[(base_name, "active")]))
        elif (base_name, "non-active") in all_groups:
            new_class = '- (%s)' % (len(group_info[
                (base_name, "non-active")]))
        else:
            new_class = "- (0)"
        new_class = "%s%s" % (new_class, "\ " * (7 - len(new_class)))
        classifiers.append(new_class)
    if org in with_tandem:
        classifiers.append("td")
    else:
        classifiers.append("-\ ")
    return classifiers

def get_db_rec(cur_id, cur_dbs):
    cur_rec = None
    for db in cur_dbs:
        try:
            cur_rec = db[cur_id]
            break
        except KeyError:
            pass
    assert cur_rec is not None, cur_id
    return cur_rec

def ids_from_group(group_name, file_base):
    """Parse out UniProt IDs from a HTML file describing the cluster group.
    """
    uniprot_ids = []
    try:
        filename = glob.glob(file_base % group_name)[0]
    except IndexError:
        return uniprot_ids
    with open(filename) as in_handle:
        soup = BeautifulSoup(in_handle)
        rows = soup.findAll("tr")
        for row in rows:
            link = row.find("a", href=re.compile("uniprot"))
            if link:
                uniprot_ids.append(link.string)
    return uniprot_ids

def set_org_layout(gv_graph, org_labels):
    """Organize the organisms to be layed out correctly in the final graph.
    """
    org_nodes = []
    pred_nodes = []
    for node in gv_graph:
        is_org_list = [1 for o in org_includes 
                if node.replace("\\n", " ").find(o) >= 0]
        if sum(is_org_list) > 0:
            node.attr['shape'] = 'record'
            node.attr['fontsize'] = '12'
            #node.attr['style'] = 'rounded'
            #node.attr['nojustify'] = 'false'
            org_nodes.append(node)
            pred_nodes.extend(gv_graph.predecessors(node))
        else:
            node.attr['shape'] = 'plaintext'
            node.attr['fontsize'] = '10'
    sub_org = gv_graph.subgraph(org_nodes)
    sub_org.graph_attr['rank'] = 'max'
    sub_pred = gv_graph.subgraph(pred_nodes)
    sub_pred.graph_attr['rank'] = 'same'
    return gv_graph

def get_available_dbs(db_dir):
    all_dbs = []
    for fname in glob.glob(os.path.join(db_dir, '*.db')):
        base, ext = os.path.splitext(fname)
        all_dbs.append(shelve.open(base))
    return all_dbs

def build_tax_graph(cur_dbs, org_includes, org_classifications):
    """Build a taxonomy graph of all organisms represented with a given domain.

    This is a relative graph in that all distances are 1 back to the top level,
    but will allow estimation of evolutionary distances based on the 
    NCBI/UniProt organism lineage trees.
    """
    org_labels = get_org_names_padded(org_includes)
    #tax_graph = networkx.DiGraph()
    tax_graph = pgv.AGraph(directed=True)
    organisms_done = []
    use_dict = collections.defaultdict(int)
    orgs_to_do = list(org_includes)
    while len(orgs_to_do) > 0:
        cur_org = orgs_to_do.pop(-1)
        for cur_db in cur_dbs:
            for db_domain in cur_db.keys():
                db_item = cur_db[db_domain]
                if (db_item["org_scientific_name"] not in organisms_done and
                        db_item["org_scientific_name"] in [cur_org]):
                    org_lineage = [node for node in db_item["org_lineage"]
                            if node in branch_includes]
                    full_lineage = org_lineage + \
                            [get_org_label(db_item["org_scientific_name"],
                                org_labels, org_classifications)]
                    #print db_item["org_lineage"] + [db_item["org_scientific_name"]]
                    for index in range(len(full_lineage) - 1):
                        use_dict[full_lineage[index]] += 1
                        tax_graph.add_edge(full_lineage[index],
                                full_lineage[index + 1])
                    organisms_done.append(db_item["org_scientific_name"])
                    try:
                        cur_org = orgs_to_do.pop(-1)
                    except IndexError:
                        break
    return tax_graph

def get_org_names_padded(org_includes):
    max_org = max([len(o) for o in org_includes])
    label_info = {}
    for org in org_includes:
        org_one, org_two = org.split()
        label = "%s %s%s" % (org_one, org_two, "\ " * (max_org - len(org)))
        label_info[org] = label
    return label_info

def get_org_label(org, org_labels, org_classifications):
    name = org_labels[org]
    classifications = org_classifications[org]
    return "{%s|%s}" % (name, "|".join(classifications))

if __name__ == "__main__":
    main(sys.argv[1])
