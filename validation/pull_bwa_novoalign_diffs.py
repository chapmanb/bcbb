#!/usr/bin/env python
"""Extract concordant variant differences between bwa and novoalign, focusing on mapping differences.

Requires:
  bedtools, pybedtools, vcflib
"""
import os
import subprocess
import sys

import pybedtools
import yaml

def main(config_file):
   with open(config_file) as in_handle:
       config = yaml.load(in_handle)
   out_dir = config["dirs"]["work"]
   if not os.path.exists(out_dir):
       os.makedirs(out_dir)
   consub_vcf = get_concordant_subset(config["calls"]["bwa"],
                                      config["calls"]["novoalign"],
                                      config["ref"], out_dir)
   if config.get("callable"):
      nocall_vcf = get_nocallable_subset(consub_vcf, config["callable"]["novoalign"])
   else:
      nocall_vcf = consub_vcf
   orig_nocall_vcf = subset_original_vcf(nocall_vcf, config["calls"]["bwa-orig"],
                                         config["ref"])
   for fname in [consub_vcf, nocall_vcf, orig_nocall_vcf]:
       with open(fname) as in_handle:
           total = sum([1 for line in in_handle if not line.startswith("#")])
       print fname, total

def subset_original_vcf(base_vcf, orig_vcf, ref_file):
    out_file = "{base}-orig.vcf".format(base=os.path.splitext(base_vcf)[0])
    if not os.path.exists(out_file):
        cmd = "vcfintersect -i {base_vcf} -r {ref_file} {orig_vcf} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def get_nocallable_subset(base_vcf, cmp_bed):
    """Retrieve subset of calls in base_vcf not in cmp_bed.
    """
    out_file = "{base}-nocallable.vcf".format(base=os.path.splitext(base_vcf)[0])
    if not os.path.exists(out_file):
        base_bt = pybedtools.BedTool(base_vcf)
        cmp_bt = pybedtools.BedTool(cmp_bed)
        base_bt.intersect(cmp_bt, v=True).saveas(out_file + ".bt")
        with open(out_file, "w") as out_handle:
            with open(base_vcf) as in_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        out_handle.write(line)
            with open(out_file + ".bt") as in_handle:
                for line in in_handle:
                    out_handle.write(line)
    return out_file

def get_concordant_subset(base_vcf, cmp_vcf, ref_file, out_dir):
    """Retrieve subset of calls in base_vcf not in cmp_vcf.
    """
    out_file = os.path.join(out_dir, "{base}-unique.vcf"
                            .format(base=os.path.splitext(os.path.basename(base_vcf))[0]))
    if not os.path.exists(out_file):
        cmd = "vcfintersect -v -i {cmp_vcf} -r {ref_file} {base_vcf} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

if __name__ == "__main__":
    main(sys.argv[1])
