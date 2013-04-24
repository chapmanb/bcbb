#/usr/bin/env Rscript

# Plot validation results across multiple samples, broken down
# by potential discordant causes. Contains hooks to allow comparisons
# between multiple alignment and variant calling methods.
#
# Usage:
#   Rscript plot_validation_results.r <grading-summary.csv>

library(ggplot2)
library(plyr)
library(stringr)
library(scales)

args <- commandArgs(trailingOnly=TRUE)
sum_file <- args[1]

# Plot faceted breakdown of comparisons for a specific variant type
plot_variant_type <- function(d, var_type, sum_file) {
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-", var_type, ".pdf")
  
  d.type <- subset(d, d$variant.type==tolower(var_type))
  d.type$sample.caller <- factor(d.type$sample.caller, levels=rev(levels(d.type$sample.caller)))
  p <- ggplot(d.type, aes(x=sample.caller, y=value, fill=sample.caller)) +
       geom_bar(stat="identity") +
       geom_text(aes(label=value), size=4, hjust=1) +
       coord_flip() +
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
       facet_wrap(~ category, ncol=3) +
       theme(axis.title.y = element_blank()) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       scale_fill_discrete(name="Calling approach", guide = guide_legend(reverse=TRUE)) +
       ylab("Count (log scale)") +
       labs(title=paste(var_type, ":", "Variant counts per approach, by classification category"))
  ggsave(out_file, p, width=14, height=10)
}

# Plot barcharts of various prep attributes
plot_prepattr <- function(d, var_type, prep_attr, sum_file) {
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-", prep_attr, "-", var_type, ".pdf")
  
  d.type <- subset(d, d$variant.type==tolower(var_type))
  d.ready <- subset(d.type, d.type[prep_attr]!="")
  d.medians <- ddply(d.ready, c("category", prep_attr), summarize, med=median(value))
  p <- ggplot(d.ready, aes_string(x=prep_attr, y="value")) +
       geom_boxplot() +
       geom_text(data=d.medians, aes(y=med, label=round(med)), size=3, vjust=1.5) +
       facet_wrap(~ category, ncol=3) +
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))
  ggsave(out_file, p, width=10, height=10)
}

# Convert method names into prettier display output
pretty_method_names <- function(d) {
  for (level in levels(d$category)) {
    new_name <- str_c(str_split(level, "-")[[1]], collapse=" ")
    new_name <- paste0(toupper(substring(new_name, 1, 1)), substring(new_name, 2, nchar(new_name)))
    levels(d$category)[levels(d$category)==level] <- new_name
  }
  d
}

d <- read.csv(sum_file, header=TRUE)
d <- transform(d, sample.caller=paste(sample, caller))
d <- pretty_method_names(d)
print(head(d))
plot_variant_type(d, "SNP", sum_file)
plot_variant_type(d, "Indel", sum_file)
plot_prepattr(d, "SNP", "caller", sum_file)
plot_prepattr(d, "Indel", "caller", sum_file)
plot_prepattr(d, "SNP", "aligner", sum_file)
plot_prepattr(d, "Indel", "aligner", sum_file)
plot_prepattr(d, "SNP", "bamprep", sum_file)
plot_prepattr(d, "Indel", "bamprep", sum_file)
