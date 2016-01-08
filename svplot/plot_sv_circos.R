#!/usr/bin/env Rscript

# Plot circos style representation of structural variants using ggbio

library(biovizBase)
library(GenomicRanges)
library(ggbio)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(grid)
library(gridExtra)

data("CRC", package="biovizBase")
ideo <- hg19sub

args <- commandArgs(trailingOnly=TRUE)
in.file <- args[1]
out.file <- str_c(str_split(in.file, "[.]")[[1]][1], ".pdf")

sv.df = read.csv(in.file, header=T, stringsAsFactors=FALSE)
sv.gr = GRanges(sv.df$chrom1, IRanges(sv.df$start1, sv.df$end1),
                svtype=sv.df$svtype, sample=sv.df$sample,
                to.gr=GRanges(sv.df$chrom2, IRanges(sv.df$start2, sv.df$end2)))
seqlevels(sv.gr) <- seqlevels(ideo)
seqlengths(sv.gr) <- seqlengths(ideo)
seqlevels(sv.gr$to.gr) <- seqlevels(ideo)
seqlengths(sv.gr$to.gr) <- seqlengths(ideo)

print(sv.gr)
svtypes <- unique(as.character(values(sv.gr)$svtype))
print(svtypes)

sv.grl <- split(sv.gr, values(sv.gr)$sample)
crc.lst <- lapply(sv.grl, function(gr.cur){
  print(unique(as.character(values(gr.cur)$sample)))
  cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
  #cols <- RColorBrewer::brewer.pal(length(svtypes), "Set2")[length(svtypes):1]
  names(cols) <- svtypes
  p <- ggbio() + circle(gr.cur, geom="link", linked.to="to.gr", aes(color=svtype)) +
                 circle(ideo, geom="text", aes(label=seqnames), vjust=0, size=4) +
                 scale_color_manual(values = cols)  +
                 labs(title = (unique(values(gr.cur)$sample)))
})

pdf(file=out.file, width=8, height=8)
bquiet = lapply(crc.lst, print)
#arrangeGrobByParsingLegend(crc.lst, legend.idx = 1, ncol = 3)
dev.off()
