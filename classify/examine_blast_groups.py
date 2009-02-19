#!/usr/bin/env python
"""Examine BLAST hits and group into clusters based on sequence similarity.
"""
from __future__ import with_statement
import sys
import math
import os
import operator
import collections
import glob
import shelve

from Bio.Blast import NCBIXML
from Bio import Fasta
import pylab

org_includes = [
        'Homo sapiens',
        'Mus musculus',
        'Rattus norvegicus',
        'Bos taurus',
        'Gallus gallus',
        'Xenopus laevis',
        'Danio rerio',
        'Drosophila melanogaster',
        'Caenorhabditis elegans',
        'Anopheles gambiae',
        'Branchiostoma floridae',
        'Tetraodon nigroviridis',
        'Nematostella vectensis'
]

def main(blast_file):
    db_dir = os.path.join(os.getcwd(), "db")
    cur_dbs = get_available_dbs(db_dir)
    length_cutoff = 0.2
    blast_clusters, all_lengths = get_blast_clusters(blast_file, length_cutoff)
    filter_clusters = filter_by_organism(blast_clusters, org_includes, cur_dbs)
    length_plot(all_lengths, blast_file)
    cluster_grouper = SimilarityClusterGrouper(2, 200, [(0.9, 10)])
    all_groups = cluster_grouper.get_final_groups(filter_clusters)
    base, ext = os.path.splitext(blast_file)
    cluster_file = base + "-bcluster%s.txt"
    for gindex, group in enumerate(all_groups):
        print '-----------'
        with open(cluster_file % gindex, "w") as out_handle:
            for gitem in group:
                db_rec = get_db_rec(gitem, cur_dbs)
                print gitem, db_rec["org_scientific_name"]
                rec = Fasta.Record()
                rec.title = gitem
                rec.sequence = db_rec["seq"]
                out_handle.write(str(rec) + "\n")

def filter_by_organism(blast_clusters, org_includes, cur_dbs):
    """Filter the listing of clusters to include only desired organisms.
    """
    filter_clusters = collections.defaultdict(lambda: [])
    for main, related_list in blast_clusters.items():
        db_rec = get_db_rec(main, cur_dbs)
        if db_rec["org_scientific_name"] in org_includes:
            for related in related_list:
                db_rec = get_db_rec(related, cur_dbs)
                if db_rec["org_scientific_name"] in org_includes:
                    filter_clusters[main].append(related)
    return dict(filter_clusters)

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

def get_available_dbs(db_dir):
    all_dbs = []
    for fname in glob.glob(os.path.join(db_dir, '*.db')):
        base, ext = os.path.splitext(fname)
        all_dbs.append(shelve.open(base))
    return all_dbs

def get_blast_clusters(blast_file, length_cutoff):
    """Retrieve a dictionary of related items from a self-vs-self BLAST.

    The length cutoff species the parameters at which to include hits
    as being included. A list of all lengths is also found for evaluation
    and threshold setting purposes.
    """
    all_lengths = []
    blast_clusters = collections.defaultdict(lambda: [])
    with open(blast_file) as blast_handle:
        records = NCBIXML.parse(blast_handle)
        for brec in records:
            for a in brec.alignments:
                cur_lengths = []
                for hsp in a.hsps:
                    try:
                        gaps = float(hsp.gaps)
                    except TypeError:
                        gaps = 0
                    hsp_length = float(len(hsp.query)) - gaps
                    cur_lengths.append(hsp_length / float(brec.query_letters))
                all_lengths.extend(cur_lengths)
                if max(cur_lengths) > length_cutoff and a.hit_def != brec.query:
                    blast_clusters[brec.query].append(a.hit_def)
    return dict(blast_clusters), all_lengths

class SimilarityClusterGrouper:
    """Group together items into a cluster using a similarity heuristic.

    Given a dictionary where:
        key - base item
        value - list of related items, ordered by relevance

    the class calculates a set of related groups. The groups are controlled
    through a maximum and minimum group size and a set of filtering parameters
    which control how we all related items through into groups.
    """
    def __init__(self, min_group_size, max_group_size, filter_params=None):
        self._min_group = min_group_size
        self._max_group = max_group_size
        if filter_params:
            self._filter_params = filter_params
        else:
            self._filter_params = [(0.08, 2), (0.10, 2), (0.11, 3),
                    (0.125, 3), (0.15, 4), (0.2, 5), (0.9, 10)]

    def get_final_groups(self, related_items):
        """Retrieve related groups for the provided relatedness dictionary.

        This works through a set of increasingly less stringent filtering
        parameters.
        """
        filter_params = self._filter_params[:]
        final_groups = []
        while len(related_items) > 0:
            if len(filter_params) == 0:
                raise ValueError("Ran out of parameters before finding groups")
            cur_thresh, cur_related = filter_params.pop(0)
            while 1:
                filt_related = self._filter_related(related_items, cur_thresh,
                        cur_related)
                groups = self._groups_from_related_dict(filt_related)
                new_groups, related_items = self._collect_new_groups(
                        related_items, groups)
                final_groups.extend(new_groups)
                if len(new_groups) == 0 or len(related_items) == 0:
                    break
            if len(related_items) < self._max_group and len(related_items) > 0:
                final_groups.append(related_items.keys())
                related_items = {}
        return final_groups

    def _collect_new_groups(self, related_items, groups):
        """Collect new groups within our parameters, updating the ones to find.
        """
        final_groups = []
        for group_items in groups:
            final_items = [i for i in group_items if related_items.has_key(i)]
            if (len(final_items) >= self._min_group and
                    len(final_items) <= self._max_group):
                final_groups.append(final_items)
                for item in final_items:
                    del related_items[item]
        final_related_dict = {}
        for item, related in related_items.items():
            final_related = [r for r in related if related_items.has_key(r)]
            final_related_dict[item] = final_related
        return final_groups, final_related_dict

    def _filter_related(self, inital_dict, overrep_thresh=0.125, related_max=3):
        """Filter a dictionary of related terms to limit over-represented items.

        We may have some items represented many times in a set, which will lead
        to a non-useful huge cluster of related items. These are filtered,
        based on overrep_thresh, and the total items for any item is limited to
        related_max.
        """
        final_dict = {}
        all_vals = reduce(operator.add, inital_dict.values())
        for item_id, item_vals in inital_dict.items():
            final_vals = [val for val in item_vals if 
                float(all_vals.count(val)) / len(inital_dict) <= overrep_thresh]
            if len(final_vals) > 0:
                final_dict[item_id] = final_vals[:related_max]
        return final_dict

    def _groups_from_related_dict(self, related_dict):
        """Create a list of unique groups from a dictionary of relations.
        """
        cur_groups = []
        all_base = related_dict.keys()
        for base_id, cur_ids in related_dict.items():
            overlap = set(cur_ids) & set(all_base)
            if len(overlap) > 0:
                new_group = set(overlap | set([base_id]))
                is_unique = True
                for exist_i, exist_group in enumerate(cur_groups):
                    if len(new_group & exist_group) > 0:
                        update_group = new_group | exist_group
                        cur_groups[exist_i] = update_group
                        is_unique = False
                        break
                if is_unique:
                    cur_groups.append(new_group)
        return [list(g) for g in cur_groups]

def length_plot(lengths, filename):
    """Prepare a histogram plot of all lengths to look for cutoffs.
    """
    base, ext = os.path.splitext(filename)
    pylab.hist(lengths, 100)
    pylab.savefig("%s-lengths.png" % base)
    #pylab.show()

if __name__ == "__main__":
    main(sys.argv[1])
