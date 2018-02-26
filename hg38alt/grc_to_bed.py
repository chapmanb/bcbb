#!/usr/bin/env python
"""Convert GRC download into ready to apply hg38 BED file.
"""
import csv
import os

import pandas as pd

in_file = "GRCh38p12.tsv"
fai_file = "hg38.fa.fai"
out_file = "%s.bed" % os.path.splitext(in_file)[0]

with open(fai_file) as in_handle:
    chroms = {}
    for chrom in (l.split()[0] for l in in_handle):
        parts = chrom.split("_")
        if len(parts) > 1:
            chroms[parts[1]] = chrom

df = pd.read_csv(in_file, sep="\t", header=0)
df.columns = [x.replace("#", "").replace(" ", "_") for x in df.columns]

with open(out_file, "w") as out_handle:
    writer = csv.writer(out_handle, delimiter="\t")
    for row in df.itertuples():
        print(row)
        seq_id = row.Sequences_IDs.split("|")[0].replace(".", "v")
        chrom = chroms.get(seq_id)
        if chrom:
            writer.writerow(["chr%s" % row.Chr, int(row.Start.replace(",", "")) - 1,
                             row.Stop.replace(",", ""), chrom])
