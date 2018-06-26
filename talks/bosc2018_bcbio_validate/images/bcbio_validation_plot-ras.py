import sys
from bcbio.variation import validateplot

title=""
validateplot.classifyplot_from_valfile(sys.argv[1], outtype="png", title=title)
