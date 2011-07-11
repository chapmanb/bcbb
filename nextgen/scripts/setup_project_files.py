#!/usr/bin/env python
"""Setup links for a project so bcbio pipeline can be used on project-specific
delivered fastq files.

Usage:
  setup_project_files.py <YAML run config> <fastq_dir> [--project_dir]

Options:

  -p, --project_dir=<project directory>   Explicitly state project directory; don't
                                          assume cwd is the correct directory
"""

import os
import sys
import glob
from optparse import OptionParser
from itertools import izip
import yaml

# XXX: Use bcbio lib ?
from scilife.pipeline.lane import get_flowcell_id
from scilife.pipeline.fastq import get_barcoded_project_files, convert_name_to_barcode_id

def main(run_info_yaml, fastq_dir, project_dir="./"):
    with open(run_info_yaml) as in_handle:
        run_info = yaml.load(in_handle)
    dirs  = dict(work_dir = project_dir)
    if os.path.exists( os.path.join(project_dir, "data", fastq_dir)):
        dirs.update(fastq_dir = os.path.join(project_dir, "data", fastq_dir))
    else:
        dirs.update(fastq_dir = fastq_dir)

    fc_name, fc_date = get_flowcell_id(run_info, dirs['fastq_dir'], glob_ext=".fastq", check_bc=False)
    config = dict(log_dir=os.path.join(project_dir, "log"), 
                  fc_name = fc_name,
                  fc_date = fc_date
                  )
    dirs.update(fc_dir = os.path.join(project_dir, "intermediate", os.path.basename(fastq_dir), 
                                                "%s_%s" %(fc_date, fc_name)))


    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        if not os.path.exists(dirs['fc_dir']):
            log.info("Creating flow cell directory %s" % (dirs['fc_dir']))
            os.makedirs(dirs['fc_dir'])
        run_main(run_info, config, dirs)

def run_main(run_info, config, dirs):
    for info in run_info:
        process_lane(info, config, dirs)

def process_lane(info, config, dirs):
    sample_name = info.get("description", "")
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)
    log.info("Processing sample: %s; lane %s; reference genome %s" %
             (sample_name, info["lane"], genome_build))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
        bc_dir = os.path.join(dirs['fc_dir'], "%s_%s_%s_barcode" % (info['lane'], config['fc_date'], config['fc_name']))
        if not os.path.exists(bc_dir):
            log.info("Making barcode directory %s" %(bc_dir))
            os.makedirs(bc_dir)
        fq = get_barcoded_project_files(multiplex, info['lane'], dirs['fastq_dir'], config['fc_name'])
        for fqpair in fq:
            fqout = convert_name_to_barcode_id(multiplex, config['fc_name'], fqpair)
            [_link_fastq_file(fq_src, fq_tgt, bc_dir) for fq_src, fq_tgt in izip(fqpair, fqout)]
    else:
        # Just write to fc_dir?
        pass
    
def _link_fastq_file(fq_src, fq_tgt, outdir):
    if fq_src is None:
        return
    tgt = os.path.join(outdir, fq_tgt)
    if os.path.exists(tgt):
        log.warn("%s already exists: link already made, not overwriting" %(tgt))
    else:
        log.info("Linking fastq file %s to %s as %s" % (os.path.basename(fq_src), outdir, fq_tgt ))
        os.symlink(fq_src, tgt)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-p", "--project_dir", dest="project_dir")
    (options, args) = parser.parse_args()
    if len(args)<1:
        print __doc__
        sys.exit()
    kwargs = dict(
        project_dir = options.project_dir if options.project_dir else os.getcwd()
        )
    main(*args, **kwargs)
