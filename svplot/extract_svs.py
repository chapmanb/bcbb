#!/usr/bin/env python
"""Extract structural variations in regions of interest into tab delimited files for circos plots.
"""
import csv
import os
import subprocess

import pysam

def main():
    samples = {"H1047R": ["H1047R-2-3", "H1047R-2-7", "H1047R-2-9", "H1047R-P"],
               "MCF7": ["MCF71_4", "MCF71_8", "MCF72_5", "MCF7P"],
               "PC9": ["PC9", "PC9_C797S_SC1", "PC9_Van", "PC9_Van_R_C797S_SCB5"]}
    chroms = {"H1047R": ["8"],
              "MCF7": ["18", "22"],
              "PC9": ["7"]}
    callers = ["lumpy", "manta"]
    depths = [10, 25]
    approaches = ["full", "subset"]
    for project, samples in samples.items():
        for depth in depths:
            for approach in approaches:
                if approach == "full":
                    cur_chroms = []
                else:
                    cur_chroms = chroms[project]
                extract_svs(project, samples, depth, cur_chroms, callers)

def extract_svs(project, samples, depth, chroms, callers):
    """Create CSV file of structural variants of interest, for Circos plots.
    """
    allowed_chroms = set([str(x) for x in range(1, 23)])
    out_dir = "sv_plots"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = os.path.join(out_dir, "%s-svs-%s-dp%s.csv" % (project, ("subset" if chroms else "full"), depth))
    region_bed = "%s-regions.bed" % os.path.splitext(out_file)[0]
    with open(out_file, "w") as out_handle:
        with open(region_bed, "w") as region_out:
            writer = csv.writer(out_handle)
            writer.writerow(["chrom1", "start1", "end1", "chrom2", "start2", "end2", "sample", "caller", "svtype"])
            for sample in samples:
                for ci, caller in enumerate(callers):
                    sv_file = os.path.join("inputs", "%s-%s.vcf.gz" % (sample, caller))
                    for p1, p2, svtype in parse_svs(sv_file, depth):
                        if len(chroms) == 0 or (p1[0] in chroms or p2[0] in chroms):
                            if p1[0] in allowed_chroms and p2[0] in allowed_chroms:
                                writer.writerow(p1 + p2 + [sample, caller, svtype])
                                if ci == 0:
                                    for (chrom, start, end) in [p1, p2]:
                                        region_out.write("%s\t%s\t%s\t%s\n" % (chrom, start - 1, end, sample))
    # filter_shared(out_file, region_bed)
    # noshared_file = "%s-noshared%s" % os.path.splitext(out_file)
    # with open(noshared_file, "w") as out_handle:
    #     writer = csv.writer(out_handle)
    #     writer.writerow(["chrom1", "start1", "end1", "chrom2", "start2", "end2", "sample", "caller", "svtype"])
    #     for row in filter_shared(out_file, region_bed):
    #         writer.writerow(row)

def filter_shared(in_file, region_bed):
    cmd = "sort -k1,1 -k2,2n {region_bed} | bedtools merge -i - -d 10000 -c 4 -o count_distinct"
    overlaps = subprocess.check_output(cmd.format(**locals()), shell=True)
    shared = 0
    noshared = 0
    for line in overlaps.split("\n"):
        if line.strip():
            if int(line.strip().split()[-1]) > 1:
                shared += 1
            else:
                noshared += 1
    print shared, noshared
    return []

def parse_svs(in_file, depth):
    bnds = {}
    min_su_bnd = depth
    min_su_other = depth
    for rec in pysam.VariantFile(in_file):
        if passes(rec):
            if rec.info["SVTYPE"] == "BND":
                if rec.info.get("SU") is None or rec.info["SU"] > min_su_bnd:
                    p = [rec.chrom.replace("chr", ""), rec.start, rec.start]
                    if rec.id in bnds:
                        p_o = bnds[rec.id]
                        yield p_o, p, rec.info["SVTYPE"]
                        del bnds[rec.id]
                    else:
                        mate_id = rec.info.get("MATEID")
                        bnds[mate_id] = p
            else:
                if rec.info.get("SU") is None or rec.info["SU"] > min_su_other:
                    p1 = [rec.chrom.replace("chr", ""), rec.start, rec.start]
                    p2 = [rec.chrom.replace("chr", ""), rec.info["END"], rec.info["END"]]
                    yield p1, p2, rec.info["SVTYPE"]

def passes(rec):
    if rec.filter.keys() == 0 or rec.filter.keys()[0] == "PASS":
        if not rec.samples[0].get("GT"):
            return True
        elif list(set(rec.samples[0].get("GT"))) != ["N"]:
            return True
    return False

if __name__ == "__main__":
    main()