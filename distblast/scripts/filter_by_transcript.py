#!/usr/bin/env python
"""Filter an output file, removing alternative transcripts based on names.

This will filter a file either based on removing alternative transcripts of the
same size or on condensing all alternative transcripts to the longest.

Usage:
    filter_by_transcript.py <org config file> <file to filter> [--condense]
"""
import sys
import os
import csv
import re
import collections
from optparse import OptionParser

import yaml
from Bio import SeqIO

def main(config_file, to_filter, condense):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    names_to_include = get_representative_txs(config['search_file'],
            condense)
    new_ext = "onetx" if condense else "nodups"
    base, ext = os.path.splitext(to_filter)
    out_file = "%s-%s%s" % (base, new_ext, ext)
    with open(to_filter) as in_handle:
        with open(out_file, "w") as out_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(reader.next())
            for parts in reader:
                if parts[0] in names_to_include:
                    writer.writerow(parts)

def get_representative_txs(in_file, condense):
    """Organize transcripts based on base names and extensions.

    Splits items (abc.1a.2) into base names (abc.1) and extensions (a.2) and
    then picks out representative transcripts to use. If condense is True,
    then the largest item is chosen. Otherwise, any duplicates of the same size
    are eliminated.
    """
    pat = re.compile("^(?P<base>\w+\.\d+)(?P<ext>([a-z]+|\.\d+)?(\.\d+)?)$")
    txs_by_size = collections.defaultdict(list)
    with open(in_file) as in_handle:
        for rec in SeqIO.parse(in_handle, "fasta"):
            match = pat.match(rec.id)
            txs_by_size[match.group("base")].append(
                    (match.group("ext"), len(rec.seq)))
    final_list = []
    for tx, choices in txs_by_size.iteritems():
        final_list.extend(_pick_representative(tx, choices, condense))
    print final_list[:10]
    return final_list

def _pick_representative(base, choices, condense):
    """Choose representative items from the list of choices based on size.
    """
    if len(choices) == 1:
        return ["".join([base, choices[0][0]])]
    else:
        cur_sizes = [-1]
        cur_choices = []
        choices.sort()
        for c, size in choices:
            # if we are condensing, we want the biggest one
            if condense:
                if size > max(cur_sizes):
                    cur_choices = [c]
            # otherwise we want non duplicates
            elif size not in cur_sizes:
                cur_choices.append(c)
            cur_sizes.append(size)
        return ["".join([base, e]) for e in cur_choices]

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-c", "--condense", dest="condense", action="store_true",
            default=False)
    options, args = parser.parse_args()
    if len(args) != 2:
        print __doc__
        sys.exit()
    main(args[0], args[1], options.condense)
