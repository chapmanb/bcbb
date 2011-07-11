#!/usr/bin/env python
"""Provide SNP and indel calling.

This script works on delivered project data. In a project 
folder sequence data is stored as fastq files in the data 
directory, grouped by flow cell like names. Optionally, 
it is possible to provide a directory containing bam 
files. 

Usage:
    exome_pipeline.py <exome pipeline YAML config file> 
                      <pruned YAML run information>
                      <project_name> <fc_dir>
                      [--project_dir=<project directory>]

Options:

  -p, --project_dir=<project directory>   Explicitly state project directory; don't
                                          assume cwd is the correct directory


The exome pipeline is similar to the post_processing.yaml file in that it holds
information about programs, algorithms and analyses. The contents are more loosely
defined than for post_processing.

The pruned YAML run information file contains run_info about the samples in this 
project only, and also contains the necessary information for conversion between 
barcode ids and barcode names, should the lanes have been multiplexed.

LANE_DATE_FC[_BCI][_1/2]_fastq.txt -> LANE_DATE_FC[_RUNINFONAME][_1/2].fastq

In general, data is delivered as fastq files to a project. Generally, we adopt 
the following convention:

j_doe_00_01/data/fastq_dir

where j_doe_00_01 is the project name to be provided at the command line. Data 
analyses are then performed in 

j_doe_00_01/intermediate/fastq_dir/

which in turn contains subdirectories for flowcells, alignments etc.

The <fc_dir> names a flowcell directory:

j_doe_00_01/intermediate/fastq_dir/fc_dir

in which the delivered fastq files have been renamed to their original names
and link back to the files in j_doe_00_01/data/fastq_dir directory. Hence,
the bcbio modules can be directly applied to the file names in fc_dir.

Relinking of files is done in the script setup_project_files.py.

Requires:
  - bcftools
  - GATK
  - annovar
"""


import os
import sys
from optparse import OptionParser

import yaml

import subprocess
import glob
import collections

from scilife.log import create_log_handler
from scilife.pipeline import log
from scilife.pipeline.lane import make_lane_items, get_flowcell_id

from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio import utils
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.broad import BroadRunner
from bcbio.ngsalign import bwa
from bcbio.pipeline import lane
from bcbio.pipeline import sample
from bcbio.pipeline.merge import organize_samples

def main(config_file, fastq_dir, run_info_yaml=None, project_dir=None):
    if project_dir == None:
        project_dir = os.getcwd()
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        run_main(config, config_file, fastq_dir, project_dir, run_info_yaml)

def run_main(config, config_file, fastq_dir, project_dir, run_info_yaml):
    work_dir = project_dir
    (_, fq_name) = os.path.split(fastq_dir)
    align_dir = os.path.join(work_dir, "intermediate",  fq_name, "alignments")

    run_info = _get_run_info(None, None, config, run_info_yaml)
    fc_name, fc_date = get_flowcell_id(run_info['details'], fastq_dir, check_bc=False,  glob_ext=".fastq")

    fastq_dir, galaxy_dir, config_dir = _get_full_paths(fastq_dir, config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = dict(fastq = fastq_dir, galaxy= galaxy_dir, 
                align = align_dir, 
                work = os.path.join(project_dir, "intermediate", os.path.basename(fastq_dir), "%s_%s" %(fc_date, fc_name)),
                config = config_dir, flowcell = None, 
                fc_dir = os.path.join(project_dir, "intermediate", os.path.basename(fastq_dir), "%s_%s" %(fc_date, fc_name))
                )
                
    # Since demultiplexing is already done, just extract run_items
    run_items = run_info['details']
    lane_items = []
    for info in run_items:
        lane_items.extend(make_lane_items(info, fc_date, fc_name, dirs, config))

    _run_parallel("process_alignment", lane_items, dirs, config)
    
    # Process samples
    sample_files, sample_fastq, sample_info = \
                  organize_samples(dirs, fc_name, fc_date, run_items)
    samples = ((n, sample_fastq[n], sample_info[n], bam_files, dirs, config, config_file)
               for n, bam_files in sample_files)
    # TODO: For some reason, in bcbio.broad.picardrun there is a line 
    # base = base.replace(".", "-") that screws up if one has a directory with a . in it
    _run_parallel("process_sample", samples, dirs, config)

def _get_run_info(fc_name, fc_date, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        return dict(details=run_details, run_id="")
    else:
        log.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        return galaxy_api.run_details(fc_name, fc_date)


def _run_parallel(fn_name, items, dirs, config):
    """Process a supplied function: single, multi-processor or distributed.
    """
    parallel = config["algorithm"]["num_cores"]
    if str(parallel).lower() == "messaging":
        runner = messaging.runner(dirs, config)
        return runner(fn_name, items)
    else:
        out = []
        fn = globals()[fn_name]
        with utils.cpmap(int(parallel)) as cpmap:
            for data in cpmap(fn, items):
                if data:
                    out.extend(data)
        return out

@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def process_sample(*args):
    return sample.process_sample(*args)

def _get_lane_info_from_fastq(fastq):
    fh = os.path.basename(fastq)
    lane = fh.split("_")[0]
    pu = "_".join(fh.split("_")[0:3])
    return lane, pu
    
def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-p", "--project_dir", dest="project_dir")
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict(
        project_dir = options.project_dir
        )
    main(*args, **kwargs)
