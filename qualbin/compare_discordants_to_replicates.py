#!/usr/bin/env python
"""Compare discordant variants between quality binning to replicates.

The goal is to identify if discordant variants are highly trusted or
potential false positives so we get a sense of how much real variation
we might be missing.

Usage:
    compare_discordants_to_replicates.py <config_yaml>
"""
import collections
import glob
import operator
import os
import sys

import yaml

def main(config):
    discordants = _get_all_discordants(config["dir"]["discordant"])
    replicates = _read_replicate_locs(config["dir"]["replicate"])
    filtered = _get_filtered(config["dir"]["orig"], config["orig"])
    all_reasons = []
    for (name1, name2), fname in discordants.iteritems():
        locs = collections.defaultdict(int)
        reasons = collections.defaultdict(int)
        other_discordant_locs = read_coords(discordants[(name2, name1)], only_pos=True)
        with open(fname) as in_handle:
            for line in (x for x in in_handle if not x.startswith("#")):
                coords = coords_from_line(line)
                cur_type = "missing"
                for rep_name, rep_locs in replicates.iteritems():
                    if coords in rep_locs:
                        cur_type = rep_name
                        break
                locs["total"] += 1
                locs[cur_type] += 1
                if cur_type.endswith("concordance"):
                    filter_name = filtered[name2].get(coords)
                    dp = get_info_item(line, "DP")
                    if filter_name is not None:
                        reasons[filter_name[0].replace("GATKStandard", "")] += 1
                    elif tuple(coords[:2]) in other_discordant_locs:
                        reasons["het/hom/indel"] += 1
                    elif dp < 10:
                        reasons["low_depth"] += 1
                    elif dp < 25:
                        reasons["mod_depth"] += 1
                    else:
                        if "allbin" in [name1, name2]:
                            print line.strip()
                        reasons["other"] += 1
        print name1, name2, dict(reasons), dict(locs)
        if "allbin" in [name1, name2]:
            all_reasons.append(dict(reasons))
    combine_reasons(all_reasons)

def combine_reasons(xs):
    final = xs[0]
    for x in xs[1:]:
        for k, v in x.iteritems():
            if not k.startswith("het/hom"):
                try:
                    final[k] += v
                except KeyError:
                    final[k] = v
    items = sorted(final.iteritems(), key=operator.itemgetter(1), reverse=True)
    for k, v in items:
        print k, v

def get_info_item(line, name):
    info_parts = line.split("\t")[7].split(";")
    item_part = [x for x in info_parts if x.startswith("%s=" % name)][0]
    _, item = item_part.split("=")
    return float(item)

def _get_all_discordants(dname):
    out = {}
    for fname in glob.glob(os.path.join(dname, "*discordance.vcf")):
        _, name1, name2, _ = os.path.basename(fname).split("-")
        if "std" in [name1, name2]:
            out[(name1, name2)] = fname
    return out

def _get_filtered(dname, orig_by_name):
    out = {}
    for name, fname in orig_by_name.iteritems():
        filtered = {}
        with open(os.path.join(dname, fname)) as in_handle:
            for line in (x for x in in_handle if not x.startswith("#")):
                if line.split("\t")[6] != "PASS":
                    filtered[coords_from_line(line)] = (line.split("\t")[6],
                                                        get_info_item(line, "MQ"),
                                                        get_info_item(line, "QD"))
        out[name] = filtered
    return out

def _read_replicate_locs(dname):
    out = {}
    for fname in glob.glob(os.path.join(dname, "*.vcf")):
        name = os.path.splitext(os.path.basename(fname))[0]
        out[name] = read_coords(fname)
    return out

def coords_from_line(line, only_pos=False):
    parts = line.split("\t")
    if only_pos:
        return tuple(parts[:2])
    else:
        return tuple(parts[:2] + parts[3:5])

def read_coords(f, only_pos=False):
    coords = []
    with open(f) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                coords.append(coords_from_line(line, only_pos))
    return set(coords)

if __name__ == "__main__":
    with open(sys.argv[1]) as in_handle:
        config = yaml.load(in_handle)
    main(config)
