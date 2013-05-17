#/usr/bin/env Rscript

# Plot validation results for comparisons of aligner parameters.
#
# Usage:
#   Rscript plot_alignerparam_comparisons.r <grading-summary-prep.csv>

library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
sum_file <- args[1]

# Convert method names into prettier display output
pretty_method_names <- function(d) {
  for (level in levels(d$category)) {
    new_name <- str_c(str_split(level, "-")[[1]], collapse=" ")
    new_name <- paste0(toupper(substring(new_name, 1, 1)), substring(new_name, 2, nchar(new_name)))
    levels(d$category)[levels(d$category)==level] <- new_name
  }
  d
}

plot_aligner_parameters <- function(d, sum_file) {
  cats.todo <- c("Concordant", "Discordant extra total", "Discordant missing total", "Discordant shared total")
  d <- d[d$category %in% cats.todo,]
  d$variant.type <- factor(d$variant.type, levels=c("snp", "indel"))
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-alignerparams.png")

  p <- ggplot(d, aes(x=sample, y=value.floor)) +
       geom_bar(stat="identity") +
       geom_text(data=d, aes(y=value.floor + 2500, label=value), size=3, vjust=1.5) +
       facet_grid(variant.type ~ category) +
       theme(axis.title.y = element_text(size=7)) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size=6, angle=15)) +
       theme(strip.text.x = element_text(size=6)) +
       theme(plot.title = element_text(size=7)) +
       ylab("Variant count") +
       labs(title="Comparison of alignment parameters with GATK UnifiedGenotyper variant calling.")
  ggsave(out_file, p, width=8, height=3)
}

d <- read.csv(sum_file, header=TRUE)
d <- pretty_method_names(d)
print(head(d))
plot_aligner_parameters(d, sum_file)
