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
