#/usr/bin/env Rscript

# Plot validation results across multiple samples, broken down
# by potential discordant causes. Contains hooks to allow comparisons
# between multiple alignment and variant calling methods.
#
# Usage:
#   Rscript plot_validation_results.r <grading-summary-prep.csv>

library(ggplot2)
library(plyr)
library(stringr)
library(scales)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
sum_file <- args[1]

# ## General plots

# Plot faceted breakdown of comparisons for a specific variant type
plot_variant_type <- function(d, var_type, sum_file) {
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-", var_type, ".pdf")
  
  d.type <- subset(d, d$variant.type==tolower(var_type))
  d.type$sample.caller <- factor(d.type$sample.caller, levels=rev(levels(d.type$sample.caller)))
  colors <- brewer.pal(8, "Dark2")
  pal <- colorRampPalette(colors)
  p <- ggplot(d.type, aes(x=sample.caller, y=value, fill=sample.caller)) +
       geom_bar(stat="identity") +
       geom_text(aes(label=value), size=4, hjust=1) +
       coord_flip() +
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
       facet_wrap(~ category, ncol=3) +
       theme(axis.title.y = element_blank()) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       scale_fill_manual(name="Calling approach", guide = guide_legend(reverse=TRUE),
                         values=pal(length(levels(d$sample.caller)))) +
       ylab("Count (log scale)") +
       labs(title=paste(var_type, ":", "Variant counts per approach, by classification category"))
  ggsave(out_file, p, width=14, height=10)
}

# Plot barcharts of various prep attributes
plot_prepattr <- function(d, var_type, prep_attr, sum_file) {
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-", prep_attr, "-", var_type, ".pdf")
  
  d.type <- subset(d, d$variant.type==tolower(var_type))
  d.ready <- subset(d.type, d.type[prep_attr]!="")
  d.medians <- ddply(d.ready, c("category", prep_attr), summarize, med=median(value), med.floor=median(value.floor) + 2000)
  p <- ggplot(d.ready, aes_string(x=prep_attr, y="value.floor")) +
       geom_boxplot() +
       geom_text(data=d.medians, aes(y=med.floor, label=round(med)), size=3, vjust=1.5) +
       facet_wrap(~ category, ncol=3) + 
       theme(axis.title.y = element_blank()) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       labs(title=paste(var_type, ":", "Variant counts, by classification category", "for", prep_attr))
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

# ## Custom plots

plot_aligner_differences <- function(d) {
  cats.todo <- c("Concordant", "Discordant extra total", "Discordant missing total", "Discordant shared total")
  d <- subset(d, d$bamprep=="gatk")
  d <- subset(d, d$caller=="gatk")
  d <- d[d$category %in% cats.todo,]
  d$variant.type <- factor(d$variant.type, levels=c("snp", "indel"))
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-alignerdiff.png")

  p <- ggplot(d, aes(x=aligner, y=value.floor)) +
       geom_bar(stat="identity") +
       geom_text(data=d, aes(y=value.floor + 2500, label=value), size=3, vjust=1.5) +
       facet_grid(variant.type ~ category) +
       theme(axis.title.y = element_text(size=7)) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size=7)) +
       theme(strip.text.x = element_text(size=6)) +
       theme(plot.title = element_text(size=7)) +
       ylab("Variant count") +
       labs(title="Alignment with GATK post-alignment preparation and GATK UnifiedGenotyper variant calling")
  ggsave(out_file, p, width=5, height=3)
}

plot_bamprep_differences <- function(d) {
  cats.todo <- c("Concordant", "Discordant extra total", "Discordant missing other", "Discordant missing total")
  d <- subset(d, d$aligner=="bwa")
  d <- subset(d, d$caller=="gatk")
  d <- d[d$category %in% cats.todo,]
  d$variant.type <- factor(d$variant.type, levels=c("snp", "indel"))
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-bamprepdiff.png")

  p <- ggplot(d, aes(x=bamprep, y=value.floor)) +
       geom_bar(stat="identity") +
       geom_text(data=d, aes(y=value.floor + 2500, label=value), size=3, vjust=1.5) +
       facet_grid(variant.type ~ category) +
       theme(axis.title.y = element_text(size=7)) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size=7)) +
       theme(strip.text.x = element_text(size=6)) +
       theme(plot.title = element_text(size=7)) +
       ylab("Variant count") +
       labs(title="Comparing post-alignment preparation with bwa alignment and GATK UnifiedGenotyper variant calling")
  ggsave(out_file, p, width=5, height=3)
}

plot_caller_differences <- function(d) {
  cats.todo <- c("Concordant", "Discordant extra low coverage", "Discordant missing low coverage", "Discordant shared hethom")
  d <- subset(d, d$aligner=="bwa")
  d <- subset(d, d$bamprep=="gatk")
  d <- d[d$category %in% cats.todo,]
  d$variant.type <- factor(d$variant.type, levels=c("snp", "indel"))
  out_file <- str_c(str_split(sum_file, "[.]")[[1]][1], "-callerdiff.png")

  p <- ggplot(d, aes(x=caller, y=value.floor)) +
       geom_bar(stat="identity") +
       geom_text(data=d, aes(y=value.floor + 2500, label=value), size=2, vjust=1.5) +
       facet_grid(variant.type ~ category) +
       theme(axis.title.y = element_text(size=7)) +
       theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size=6, angle=20)) +
       theme(strip.text.x = element_text(size=6)) +
       theme(plot.title = element_text(size=8)) +
       ylab("Variant count") +
       labs(title="Comparing calling methods with bwa alignment and GATK post-alignment preparation")
  ggsave(out_file, p, width=6, height=3)
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

plot_aligner_differences(d)
plot_bamprep_differences(d)
plot_caller_differences(d)
