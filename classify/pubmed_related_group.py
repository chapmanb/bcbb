#!/usr/bin/env python
"""Group a list of PubMed IDs based on their relationships in Entrez.

This is a demonstration script tackling the problem of organizing a list
of journal articles into most related groups. The input is a list of pubmed IDs
and the resulting groups are printed, demonstrating the functionality.

It works by trying to group together most closely related items first and
gradually relaxing stringency until all items are placed into related groups.

Usage:
    pubmed_related_group.py <pmids>
"""
from __future__ import with_statement
import sys
import operator

from Bio import Entrez

def main(pmid_file):
    pmids = []
    with open(pmid_file) as in_handle:
        for line in in_handle:
            pmids.append(line.strip())
    entrez_grouper = EntrezRelatedGrouper(3, 10)
    all_groups = entrez_grouper.get_pmid_groups(pmids)
    print all_groups

class EntrezRelatedGrouper:
    """Group journal articles using the Entrez Elink related query.
    """
    def __init__(self, min_group_size, max_group_size):
        self._min_group = min_group_size
        self._max_group = max_group_size

        self._filter_params = [(0.08, 2), (0.10, 2), (0.11, 3),
                (0.125, 3), (0.15, 4), (0.2, 5), (0.9, 10)]

    def get_pmid_groups(self, pmids):
        """Retrieve related groups for the passed article PubMed IDs.

        This works through a set of increasingly less stringent filtering
        parameters, placing all PubMed IDs in groups based on related articles
        from Entrez.
        """
        pmid_related = self._get_elink_related_ids(pmids)
        filter_params = self._filter_params[:]
        final_groups = []
        while len(pmid_related) > 0:
            if len(filter_params) == 0:
                raise ValueError("Ran out of parameters before finding groups")
            cur_thresh, cur_related = filter_params.pop(0)
            while 1:
                filt_related = self._filter_related(pmid_related, cur_thresh,
                        cur_related)
                groups = self._groups_from_related_dict(filt_related)
                new_groups, pmid_related = self._collect_new_groups(
                        pmid_related, groups)
                final_groups.extend(new_groups)
                if len(new_groups) == 0:
                    break
            if len(pmid_related) < self._max_group:
                final_groups.append(pmid_related.keys())
                pmid_related = {}
        return final_groups

    def _collect_new_groups(self, pmid_related, groups):
        """Collect new groups within our parameters, updating the ones to find.
        """
        final_groups = []
        for group_items in groups:
            final_items = [i for i in group_items if pmid_related.has_key(i)]
            if (len(final_items) >= self._min_group and
                    len(final_items) <= self._max_group):
                final_groups.append(final_items)
                for item in final_items:
                    del pmid_related[item]
        final_related_dict = {}
        for pmid, related in pmid_related.items():
            final_related = [r for r in related if pmid_related.has_key(r)]
            final_related_dict[pmid] = final_related
        return final_groups, final_related_dict

    def _get_elink_related_ids(self, pmids):
        """Query Entrez elink for pub med ids related to the passed list.

        Returns a dictionary where the keys are input pubmed ids and the keys
        are related PubMed IDs, sorted by score.
        """
        pmid_related = {}
        for pmid in pmids:
            handle = Entrez.elink(dbform='pubmed', db='pubmed', id=pmid)
            record = Entrez.read(handle)
            cur_ids = []
            for link_dict in record[0]['LinkSetDb'][0]['Link']:
                cur_ids.append((int(link_dict.get('Score', 0)),
                    link_dict['Id']))
            cur_ids.sort()
            cur_ids.reverse()
            local_ids = [x[1] for x in cur_ids if x[1] in pmids]
            if pmid in local_ids:
                local_ids.remove(pmid)
            pmid_related[pmid] = local_ids
        return pmid_related

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

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print __doc__
        sys.exit()
    main(sys.argv[1])
