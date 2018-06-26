import sys
from bcbio.variation import validateplot

title="DeepVariant: CHM haploid diploid"
validateplot.classifyplot_from_valfile(sys.argv[1], outtype="png", title=title)
