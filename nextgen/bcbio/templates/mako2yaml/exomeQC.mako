## Configuration for a single-end run
## Analysis files
target: ${target}
baits: ${bait}
bedfile: ${bedfile}

## Sample label for easy reference
label: ${label}

## Genome: one of NA, hg18, hg19
genome: ${genome}

## Save configuration
outdir: ${outdir}
pairedend: ${pairedend}
saverda: TRUE

## Analyses - either TRUE or FALSE
analyses:
  barplot: TRUE
  specificity: TRUE
  coverage: TRUE
  enrichment: TRUE
  duplicates: TRUE
