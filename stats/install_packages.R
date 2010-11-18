# Automate installation of useful packages when upgrading R.
# This contains everything we normally want with R so we can
# make updates smoother.
# 
# Run with:
#   Rscript install_packages.R

# Set mirrors
cran.repos <- getOption("repos")
cran.repos["CRAN" ] <- "http://software.rc.fas.harvard.edu/mirrors/R/"
options(repos=cran.repos)
source("http://bioconductor.org/biocLite.R")


repo.installer <- function(repos, install.fn) {
  update.or.install <- function(pname) {
    if (pname %in% installed.packages())
      update.packages(lib.loc=c(pname), repos=repos, ask=FALSE)
    else
      install.fn(pname)
  }
}

# standard packages
std.pkgs <- c("ggplot2", "rjson", "sqldf", "NMF", "caTools", "ape", "snowfall",
	      "multicore", "mclust")
std.installer = repo.installer(cran.repos, install.packages)
lapply(std.pkgs, std.installer)
# bioconductor packages
bioc.pkgs <- c("Biostrings", "ShortRead", "BSgenome", "edgeR", "GEOquery",
	       "GOstats", "rtracklayer", "biomaRt", "Rsamtools", "PICS", "MotIV", 
	       "rGADEM", "GO.db", "KEGG.db", "org.Hs.eg.db", "org.Mm.eg.db",
	       "affy", "affydata", "affyio", "celeganscdf", "hgu95av2cdf",
	       "preprocessCore", "ChIPpeakAnno")
bioc.installer = repo.installer(biocinstallRepos(), biocLite)
lapply(bioc.pkgs, bioc.installer)

# update anything which was not explicitly specified
update.packages(repos=biocinstallRepos(), ask=FALSE)
update.packages(ask=FALSE)
