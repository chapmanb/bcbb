## Plot alignment changes resulting from Illimina quality binning
## Usage:
##   Rscript plot_alignment_qualbin.R orig_bam qualbinned_bam

library(ggplot2)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
orig_bam <- args[1]
qualbin_bam <- args[2]
plot_file <- "qualbin-alignment-changes.png"

get_idxstats <- function(in_file) {
  stats <- system(paste("samtools idxstats", in_file), intern=TRUE)
  stats_df <- read.table(text=stats, sep="\t", col.names=c("contig", "length", "mapped", "unmapped"))
}

orig_idxstats <- get_idxstats(orig_bam)
qualbin_idxstats <- get_idxstats(qualbin_bam)

idxstats <- merge(orig_idxstats, qualbin_idxstats, by="contig")
total <- sum(idxstats$mapped.x)
idxstats$change <- abs(ifelse(idxstats$mapped.x > 0,
                              idxstats$mapped.y - idxstats$mapped.x,
                              idxstats$unmapped.y - idxstats$unmapped.x))
print(head(idxstats))
idxstats <- subset(idxstats, select=c("contig", "change"))
idxstats <- head(idxstats, 23)
idxstats$type <- ifelse(idxstats$contig == "*", "unmapped", "mapped")
idxstats$contig <- reorder(idxstats$contig, as.integer(idxstats$contig))
idxstats <- idxstats[order(as.integer(idxstats$contig)),]
print(idxstats)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot(idxstats, aes(x=contig, y=change, fill=type)) + geom_bar(stat="identity") +
     theme(axis.title.x = element_blank()) +
     theme(axis.title.y = element_blank()) +
     labs(title=paste("Read changes,",
                      as.integer(round(total / 1e6)),
                      "million total reads")) +
     scale_fill_manual(values=cbPalette)
  
ggsave(plot_file, p, width=8, height=5)
