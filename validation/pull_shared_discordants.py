#!/usr/bin/env python
"""Extract discordant variants found in multiple calling methods.

These are potential incorrect calls in the reference materials.

Requires:
  vcfintersect
"""
import itertools
import os
import subprocess
import sys

import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    out_dir = config["dirs"]["work"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    call_cmps = get_call_cmps(config["calls"], config["ref"], out_dir)
    union_file = get_union(call_cmps, config["ref"], out_dir)
    for f in list(config["calls"].itervalues()) + call_cmps + [union_file]:
       print f, variant_count(f)

def variant_count(fname):
    with open(fname) as in_handle:
        return sum([1 for line in in_handle if not line.startswith("#")])

def get_union(call_cmps, ref_file, out_dir):
   cur_union = call_cmps[0]
   for i, next_cmp in enumerate(call_cmps[1:]):
       out_file = os.path.join(out_dir, "union-{i}.vcf".format(i=i))
       cur_union = combine_two_vcfs(cur_union, next_cmp, ref_file, out_file)
   return cur_union

def combine_two_vcfs(f1, f2, ref_file, out_file):
    if not os.path.exists(out_file):
        cmd = "vcfintersect -u {f1} -r {ref_file} {f2} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def get_call_cmps(calls, ref_file, out_dir):
   """Retrieve pairwise intersections between all pairs of variants.
   """
   return [intersect_two_vcfs(c1, c2, calls[c1], calls[c2], ref_file, out_dir)
           for c1, c2 in itertools.combinations(calls.keys(), 2)]

def intersect_two_vcfs(n1, n2, f1, f2, ref_file, out_dir):
    out_file = os.path.join(out_dir, "{n1}-{n2}-intersect.vcf".format(**locals()))
    if not os.path.exists(out_file):
        cmd = "vcfintersect -i {f1} -r {ref_file} {f2} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

if __name__ == "__main__":
    main(sys.argv[1])
