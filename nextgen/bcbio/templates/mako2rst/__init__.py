"""
Mako templates for rst output
"""

import sys
from texttable import *

def program_info(proj_conf):
    d = proj_conf['program']
    tab = Texttable()
    tab.set_cols_align(["l", "l"])
    tab.set_cols_align(["b", "b"])
    tab.header(["Program", "Value"])
    for k in d.keys():
        tab.add_row([k, d[k]])
    return tab.draw()

# Has 19 columns, need better way of representing it
def align_metrics(d):
    if d==None:
        return
    tab = Texttable()
    first = True
    for fc in d.keys():
        md = d[fc]
        for lab in md.keys():
            data = md[lab]
            for row in data:
                if row[0] == "CATEGORY": 
                    if first:
                        tab.header(row)
                        tw = []
                        for td in row:
                            tw.append(len(td))
                        tab.set_cols_width(tw)
                        first = False
                else:
                    tab.add_row(row)
    return tab.draw()

def teqc_graphics(d, which="chrom-barplot", width="45%"):
    res = []
    if d==None:
        return
    for fc in d.keys():
        flowcell = d[fc]
        png = flowcell[which]
        if len(png) == 0:
            next
        else:
            print png
            for grf in png:
                res.append(".. figure:: %s\n    :width: %s\n\n" % (grf, width))
    return "\n".join(res)
