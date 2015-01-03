#!/bin/bash
# Provide BED files of GRCh37 regions that are problematic for structural variant calling
set -e

# Create BED file of no mappability regions to exclude from structural variant
# detection.
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/
# http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
bigWigToBedGraph wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.bedg
 cat wgEncodeCrgMapabilityAlign100mer.bedg \
     | awk -v maxmap=0.05 '{ if ($4 <= maxmap) {print $0"\tnomap"} else {print $0"\tmap"}}' \
     | bedtools groupby -g 1,5 -c 1,2,3,5 -o first,first,last,first \
     | cut -f 3-6 \
     | awk '{if ($4 == "nomap" && $3 - $2 > 150) print}' \
     | sed 's/^chr//' \
     > wgEncodeCrgMapabilityAlign100mer-nomap.bed

# Create BED file of low/high GC regions
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/gc5Base/
# http://wiki.bits.vib.be/index.php/Create_a_GC_content_track

# cut -f 1,2 /human/GRCh37/seq/GRCh37.fa.fai > GRCh37.sizes
# bedtools makewindows -g GRCh37.sizes -w 50 > GRCh37_windows50bp.bed
# bedtools nuc -fi /human/GRCh37/seq/GRCh37.fa -bed GRCh37_windows50bp.bed | cut -f 1-3,5 | grep -v ^# \
#     > GRCh37_gcs.bed
#  cat GRCh37_gcs.bed \
#      | awk '{if ($4 < 0.35 || $4 > 0.70) {print $0"\tbadgc"} else {print $0"\tokgc"}}' \
#      | bedtools groupby -g 1,5 -c 1,2,3,5 -o first,first,last,first \
#      | cut -f 3-6 \
#      | awk '{if ($4 == "badgc" && $3 - $2 >= 150) print}' \
#      > GRCh37_gcs-problem.bed

# Combine with low complexity regions into one global exclusion file

cat <(zcat LCR.bed.gz | awk '{print $0"\tlcr"}') \
    wgEncodeCrgMapabilityAlign100mer-nomap.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o distinct > sv_exclude.bed
