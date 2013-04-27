"""Prepare subset regions of full NIST NA12878 reference materials for evaluation.

Allows preparation of exome or targeted reference materials from
the full NIST NA12878 genome.

Requires:
  vcflib: https://github.com/ekg/vcflib
  bedtools: http://bedtools.readthedocs.org/en/latest/

Usage:
  prep_rm_subset.py <input_config.yaml>
"""
import os
import sys
import subprocess

import yaml
import pybedtools

def main(config_file):
    config = load_config(config_file)
    config["out_base"] = os.path.join(config["dirs"]["rm"],
                                      config["subset"]["name"])
    region_bed = intersect_beds(config["subset"]["interval"],
                                config["rm"]["interval"], config)
    final_vcf = combine_subset_vcfs(config["rm"]["vcfs"],
                                    config["rm"]["ref"],
                                    region_bed, config)
    filter_vcf(final_vcf)

def filter_vcf(in_vcf):
    out_vcf = "%s-pass%s" % os.path.splitext(in_vcf)
    with open(in_vcf) as in_handle:
        with open(out_vcf, "w") as out_handle:
            for line in in_handle:
                passes = False
                if line.startswith("#"):
                    passes = True
                else:
                    parts = line.split("\t")
                    if parts[6] in [".", "PASS"]:
                        passes = True
                if passes:
                    out_handle.write(line)

def combine_subset_vcfs(vcfs, ref_file, region_bed, config):
    out_file = os.path.join(config["dirs"]["rm"],
                            "%s.vcf" % config["subset"]["name"])
    tmp_files = []
    for i, vcf in enumerate(vcfs):
        tmp_out_file = "%s-%s.vcf" % (os.path.splitext(out_file)[0], i)
        cmd = "vcfintersect -b {region_bed} {vcf} > {tmp_out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
        tmp_files.append(tmp_out_file)
    # Need to generalize for multiple VCFs
    one_vcf, two_vcf = tmp_files
    cmd = "vcfintersect -r {ref_file} -u {two_vcf} {one_vcf} > {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    for tmp_file in tmp_files:
        os.remove(tmp_file)
    return out_file

def intersect_beds(base_bed, rm_bed, config):
    out_file = os.path.join(config["dirs"]["rm"],
                            "%s-regions.bed" % config["subset"]["name"])
    if not os.path.exists(out_file):
        base_bt = pybedtools.BedTool(base_bed)
        base_bt.intersect(rm_bed).saveas(out_file)
    return out_file

def load_config(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    dirs = config["dirs"]
    config["rm"]["vcfs"] = [os.path.join(dirs["rm"], x) for x in config["rm"]["vcfs"]]
    config["rm"]["interval"] = os.path.join(dirs["rm"], config["rm"]["interval"])
    config["subset"]["interval"] = os.path.join(dirs["rm"], config["subset"]["interval"])
    config["rm"]["ref"] = os.path.join(dirs["genome"], config["rm"]["ref"])
    return config

if __name__ == "__main__":
    main(sys.argv[1])
