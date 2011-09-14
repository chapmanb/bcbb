"""Functions for dealing with run_info configurations"""

import os

from bcbio.pipeline import log

def prune_run_info_by_description(run_info, desc):
    """Prune a run_info file by lane description"""
    run_info_ret = list()
    for info in run_info:
        if info['description'].find(desc) > -1:
            log.info("Found %s in run_info for lane %s, description: %s" %(desc, info['lane'], info['description'] ))
            run_info_ret.append(info)
    return run_info_ret
        
