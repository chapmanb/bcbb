#!/usr/bin/env python
"""Deliver samples to a project

Usage:
    sample_delivery.py <YAML run config> <yaml project description>
                       <flowcell directory> <project output directory>
                       [--flowcell_alias=<flowcell alias> --only_install_run_info
                       --make_delivery_report --dry_run]

The YAML run config stores information about the samples loaded on each lane. 
Given a project description and a flowcell directory, the script will look 
for all yaml lane description fields that (partially) match this description 
and installs the corresponding fastq files into the project output directory, 
creating a data subdirectory and a flowcell directory in the data directory. 

For a multiproject run_info file, only a subset of the lanes are actually used. 
The run_info file is therefore pruned, and the pruned file is output to the 
sample delivery directory.

Options:

  -a, --flowcell_alias=<flowcell alias>  By default, samples are delivered to 
                                         project_output_directory/data/FLOWCELL_ID
                                         This option changes FLOWCELL_ID
  -I, --only_install_run_info            Don't deliver samples, only install pruned
                                         run_info file
  -n, --dry_run                          Don't deliver samples, just list what will
                                         happen
  -R, --make_delivery_report             Make a delivery report
"""

import os
import sys
from optparse import OptionParser

import yaml
import glob
import shutil
from itertools import izip

from scilife.log import create_log_handler
from scilife.pipeline import log
from scilife.pipeline.run_info import prune_run_info_by_description
from scilife.pipeline.lane import get_flowcell_id
from scilife.pipeline.fastq import get_single_fastq_files, get_barcoded_fastq_files, convert_barcode_id_to_name

from bcbio import utils

def main(run_info_yaml, yaml_project_desc, fc_dir, project_outdir, 
         only_run_info, report, fc_alias=None):
    with open(run_info_yaml) as in_handle:
        run_info = yaml.load(in_handle)
    config = dict(log_dir=os.path.join(project_outdir, "log"),
                  yaml_project_desc = yaml_project_desc
                  )
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        run_info = prune_run_info_by_description(run_info, yaml_project_desc)

    dirs = dict(fc_dir=fc_dir, project_dir=project_outdir)
    fc_name, fc_date = get_flowcell_id(run_info, dirs['fc_dir'])
    config.update( fc_name = fc_name, fc_date = fc_date)
    config.update( fc_alias = "%s_%s" % (fc_date, fc_name) if not fc_alias else fc_alias)
    dirs.update(fc_delivery_dir = os.path.join(dirs['project_dir'], "data", config['fc_alias'] ))

    with log_handler.applicationbound():
        config = _make_delivery_directory(dirs, config)
        _save_run_info(run_info, dirs['fc_delivery_dir'], run_exit=only_run_info)
        run_main(run_info, config, dirs, report)

def run_main(run_info, config, dirs, report):
    for info in run_info:
        process_lane(info, dirs, config)

def process_lane(info, dirs, config):
    """Models bcbio process lane"""
    sample_name = info.get("description", "")
    genome_build = info.get("genome_build", None)
    multiplex = info.get('multiplex', None)
    log.info("Processing sample: %s; lane %s; reference genome %s" %
             (sample_name, info["lane"], genome_build))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
        fq = get_barcoded_fastq_files(multiplex, info['lane'], dirs['fc_dir'], config['fc_name'], config['fc_date'])
        for fqpair in fq:
            fqout = convert_barcode_id_to_name(multiplex, config['fc_name'], fqpair)
            [_deliver_fastq_file(fq_src, fq_tgt, config['fc_delivery_dir']) for fq_src, fq_tgt in izip(fqpair, fqout)]
    else:
        pass
        #fq1, fq2 = get_single_fastq_files(info['lane'], dirs['fc_dir'], config['fc_name'])

        
# TODO: Check for data in path so one could pass j_doe_00_01 or
# j_doe_00_01/data/flowcell_alias as project_output_dir?
def _make_delivery_directory(dirs, config):
    """Make the output directory"""
    outdir = dirs['fc_delivery_dir']
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        log.info("Creating data delivery directory %s" % (outdir))
    else:
        log.warn("%s already exists: not creating new directory" % (outdir))
    config.update(fc_delivery_dir=outdir)
    return config

def _deliver_fastq_file(fq_src, fq_tgt, outdir):
    if fq_src is None:
        return
    tgt = os.path.join(outdir, fq_tgt)
    if os.path.exists(tgt):
        log.warn("%s already exists: delivery already made, not overwriting" %(tgt))
    else:
        if not options.dry_run:
            log.info("Delivering fastq file %s to %s as %s" % (os.path.basename(fq_src), outdir, fq_tgt ))
            shutil.copyfile(fq_src, tgt)
    if options.dry_run:
        print "DRY_RUN: Delivering fastq file %s to %s as %s" % (os.path.basename(fq_src), outdir, fq_tgt)

def _save_run_info(run_info, outdir, run_exit=False):
    outfile = os.path.join(outdir, "project_run_info.yaml")
    if not options.dry_run:
        with open(outfile, "w") as out_handle:
            yaml.dump(run_info, stream=out_handle)
    else:
        print "DRY_RUN:"
        yaml.dump(run_info, stream=sys.stdout)
    if run_exit:
        sys.exit()

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-a", "--flowcell_alias", dest="fc_alias")
    parser.add_option("-I", "--only_install_run_info", dest="only_run_info", action="store_true",
                      default=False)
    parser.add_option("-R", "--make_delivery_report", dest="report", action="store_true",
                      default=False)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 4:
        print __doc__
        sys.exit()
    kwargs = dict(
        fc_alias = options.fc_alias,
        only_run_info = options.only_run_info,
        report   = options.report
        )
    main(*args, **kwargs)
