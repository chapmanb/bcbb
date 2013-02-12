## Plot causes of variant discordance due to quality binning
##
## Usage:
##  plot_variant_qualbin.R <reason CSV>

library(ggplot2)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
reason_csv <- args[1]
plot_file <- "qualbin-variant-changes.png"

reasons <- read.csv(reason_csv, header=FALSE, col.names=c("reason", "count"))
reasons$reason <- factor(reasons$reason, levels=unique(reasons$reason))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot(reasons, aes(x=reason, y=count)) + geom_bar(stat="identity") +
     xlim(rev(levels(reasons$reason))) +
     coord_flip() +
     theme(axis.title.x = element_blank()) +
     theme(axis.title.y = element_blank()) +
     labs(title="Full vs binned quality score discordant variants: causes")
     scale_fill_manual(values=cbPalette)
  
ggsave(plot_file, p, width=7, height=5)
