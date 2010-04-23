# Automate installation of useful packages when upgrading R.
# This contains everything we normally want with R so we can
# make updates smoother.
# 
# Run with:
#   Rscript install_packages.R

# Set mirror
r <- getOption("repos")
r["CRAN" ] <- "http://software.rc.fas.harvard.edu/mirrors/R/"
options(repos=r)

# standard packages
install.packages("ggplot2")
install.packages("rjson")
install.packages("sqldf")
install.packages("NMF")
install.packages("caTools")

# bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("edgeR")
biocLite("GEOquery")
biocLite("GOstats")
biocLite("rtracklayer")
biocLite("biomaRt")
biocLite("GO.db")
biocLite("KEGG.db")
biocLite("org.Hs.eg.db")
biocLite("org.Mm.eg.db")
