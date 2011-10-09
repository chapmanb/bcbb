#!/usr/bin/env python
"""Deliver samples to a project

Usage:
    sample_delivery.py <flowcell_id> <project_id>
                       [--analysis_base_dir=<analysis_base_dir>
                        --archive_base_dir=<archive_base_dir>
                        --project_base_dir=<project_base_dir
                        --flowcell_alias=<flowcell alias>
                        --project_desc=<project_desc>
                        --lanes=<lanes>
                        --copy_data
                        --only_install_run_info
                        --only_install_fastq
                        --make_delivery_report --dry_run --verbose]


Given a flowcell_id, a project_desc or a comma-separated list of lane numbers,
and a project_id, all files in <analysis_base_dir>/<flowcell_id> related to the
project defined in <archive_base_dir>/run_info.yaml will be moved to
<project_base_dir>/<project_id>.

For a multiproject run_info file, only a subset of the lanes can be used.
The run_info file is therefore pruned, and the pruned file is output to the 
sample delivery directory. The pruning is based on the options <project yaml id>,
or <lanes>. Keyword ALL delivers all lanes.

Options:

  -X, --analysis_base_dir=<analysis_base_dir>   Base dir where automated_initial_analysis.py has
                                                been run
  -A, --archive_base_dir=<archive_base_dir>     Archiving base directory
  -P, --project_base_dir=<project_base_dir>     Project base directory. Files will be installed in
                                                <project_base_dir>/<project_id>
  
  -a, --flowcell_alias=<flowcell alias>         By default, samples are moved to a directory named
                                                <flowcell_id>. This option changes output directory to
                                                <flowcell_alias>.

  -y, --project_desc=<project_desc>       Project description in description field of run_info file, or ALL.
  -l, --lanes=<lanes>                           Comma-separated list of integers corresponding to lanes

  -i, --only_install_run_info                   Only install pruned run_info file.
  -f, --only_install_fastq                      Only install fastq files.
  -r, --run_info_name                           Change run_info_name (default: run_info.yaml)
  -c, --copy_data                               Don't move data, just copy
  -n, --dry_run                                 Don't do anything samples, just list what will happen
  -v, --verbose                                 Print some more information
"""

import os
import sys
from optparse import OptionParser

import yaml
import glob
import shutil
import pprint
from itertools import izip

from bcbio.log import create_log_handler
from bcbio.pipeline import log
from bcbio.pipeline.run_info import prune_run_info_by_description
from bcbio.pipeline.lane import get_flowcell_id
from bcbio.pipeline.fastq import get_single_fastq_files, get_barcoded_fastq_files, convert_barcode_id_to_name, get_fastq_files

from bcbio import utils

pp = pprint.PrettyPrinter(indent=4)

def main(flowcell_id, project_id, fc_alias=None, project_desc=None, lanes=None):
    if project_desc is None and lanes is None:
        log.error("No project description or lanes provided: cannot deliver files without this information")
        sys.exit()
    run_info_yaml = os.path.join(options.archive_base_dir, flowcell_id, options.run_info_name)
    with open(run_info_yaml) as in_handle:
        run_info = yaml.load(in_handle)
    fc_dir = os.path.join(options.analysis_base_dir, flowcell_id)
    project_outdir = os.path.join(options.project_base_dir, project_id)

    ## Config file for sample_delivery.py
    config = dict(
        log_dir=os.path.join(project_outdir, "log"),
        project_desc = project_desc,
        lanes = lanes
        )
    
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        run_info = prune_run_info_by_description(run_info, project_desc, lanes)
    if len(run_info) == 0:
        log.error("No lanes found with matching description %s: please check your run_info.yaml file" % project_desc)
        sys.exit()

    dirs = dict(fc_dir=fc_dir, project_dir=project_outdir,
                analysis_base_dir = options.analysis_base_dir,
                archive_base_dir = options.archive_base_dir,
                project_base_dir = options.project_base_dir
                )
    fc_name, fc_date = get_flowcell_id(run_info, dirs['fc_dir'])
    config.update(fc_name = fc_name, fc_date = fc_date)
    config.update(fc_alias = "%s_%s" % (fc_date, fc_name) if not fc_alias else fc_alias)
    dirs.update(fc_delivery_dir = os.path.join(dirs['project_dir'], "data", "nobackup", config['fc_alias'] ))
    dirs.update(data_delivery_dir = os.path.join(dirs['project_dir'], "intermediate", "nobackup", config['fc_alias'] ))

    if options.verbose:
        print "=" * 50 + "\nConfiguration:\n" + "=" * 50
        pp.pprint(config)
        print "=" * 50 + "\nDirectories:\n" + "=" * 50
        pp.pprint(dirs)
        print "=" * 50 + "\nSamples to deliver:\n" + "=" * 50
        pp.pprint(run_info)

    with log_handler.applicationbound():
        config = _make_delivery_directory(dirs, config)
        _save_run_info(run_info, dirs['fc_delivery_dir'], run_exit=options.only_run_info)
        run_main(run_info, config, dirs)

def run_main(run_info, config, dirs):
    for info in run_info:
        process_lane(info, config, dirs)

def process_lane(info, config, dirs):
    """Models bcbio process lane"""
    sample_name = info.get("description", "")
    genome_build = info.get("genome_build", None)
    multiplex = info.get('multiplex', None)
    fc_link_dir = None
    log.info("Processing sample: %s; lane %s; reference genome %s" %
             (sample_name, info["lane"], genome_build))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
    fq = get_barcoded_fastq_files(multiplex, info, dirs['fc_dir'], config['fc_name'], config['fc_date'])

    ## Move data along with fastq files
    if not options.only_fastq:
        if multiplex:
            fc_link_dir = os.path.join(config['data_delivery_dir'], "%s_%s_%s_barcode" % (info['lane'], config['fc_date'], config['fc_name']))
        else:
            fc_link_dir = config['data_delivery_dir']
        _make_dir(fc_link_dir, "fastq.txt to fastq link directory")
        data, fastqc = _get_analysis_results(config, dirs, info['lane'])
        _deliver_data(data, fastqc, config['data_delivery_dir'])

    for fqpair in fq:
        fqout = convert_barcode_id_to_name(multiplex, config['fc_name'], fqpair)
        [_deliver_fastq_file(fq_src, fq_tgt, config['fc_delivery_dir'], fc_link_dir) for fq_src, fq_tgt in izip(fqpair, fqout)]

def _make_delivery_directory(dirs, config):
    """Make the output directory"""
    _make_dir(dirs['fc_delivery_dir'], "flowcell delivery")
    _make_dir(dirs['data_delivery_dir'], "data delivery")
    config.update(fc_delivery_dir=dirs['fc_delivery_dir'])
    config.update(data_delivery_dir=dirs['data_delivery_dir'])
    return config

def _make_dir(dir, label):
    if not os.path.exists(dir):
        os.makedirs(dir)
        log.info("Creating %s directory %s" % (label, dir))
    else:
        log.warn("%s already exists: not creating new directory" % (dir))

def _copy_or_move_data(src, tgt):
    if src is None:
        return
    if os.path.exists(tgt):
        log.warn("%s already exists: not overwriting" %(tgt))
        return
    if not options.dry_run:
        if options.copy:
            log.info("Copying file %s to %s" % (src, tgt))
            shutil.copyfile(src, tgt)
        else:
            log.info("Moving file %s to %s" % (src, tgt))
            shutil.move(src, tgt)
    if options.dry_run:
        if options.copy:
            print "DRY_RUN: Copying file %s to %s" % (src, tgt)
        else:
            print "DRY_RUN: Moving file %s to %s" % (src, tgt)

def _deliver_fastq_file(fq_src, fq_tgt, outdir, fc_link_dir=None):
    _copy_or_move_data(fq_src, os.path.join(outdir, fq_tgt))
    if not fc_link_dir is None:
        _link_data(os.path.join(outdir, fq_tgt), os.path.join(fc_link_dir, os.path.basename(fq_src)))

def _link_data(src, tgt):
    if src is None:
        return
    
    if options.dry_run:
        print "DRY_RUN: Linking file %s to %s" % (src, tgt)
    else:
        log.info("Linking file %s to %s" % (src, tgt))
        os.symlink(src, tgt)

def _deliver_data(data, fastqc, outdir):
    for src in data:
        tgt = os.path.join(outdir, os.path.basename(src))
        _copy_or_move_data(src, tgt)

    for src in fastqc:
        tgt = os.path.join(outdir, "fastqc", os.path.basename(src))
        _copy_or_move_data(src, tgt)
        
def _get_analysis_results(config, dirs, lane):
    """Get analysis results

    For now just glob the analysis directory for fastqc output and files with the give flowcell name
    """
    flowcell = "_".join([str(lane), config['fc_date'], config['fc_name']])
    glob_str = os.path.join(dirs['fc_dir'], flowcell + "*.*")
    data = glob.glob(glob_str)
    glob_str = os.path.join(dirs['fc_dir'], "fastqc", flowcell + "*")
    fastqc = glob.glob(glob_str)
    return data, fastqc

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
    usage = """
    sample_delivery.py <flowcell_id> <project_id>
                       [--analysis_base_dir=<analysis_base_dir>
                        --archive_base_dir=<archive_base_dir>
                        --project_base_dir=<project_base_dir
                        --flowcell_alias=<flowcell alias>
                        --project_desc=<project yaml id>
                        --lanes=<lanes>
                        --copy_data
                        --only_install_run_info
                        --only_install_fastq
                        --make_delivery_report --dry_run --verbose]

    For more extensive help type sample_delivery.py

"""

    parser = OptionParser(usage=usage)
    parser.add_option("-X", "--analysis_base_dir", dest="analysis_base_dir",
                      default="/bubo/proj/a2010002/nobackup/romanvg")
    parser.add_option("-A", "--archive_base_dir", dest="archive_base_dir",
                      default="/bubo/proj/a2010002/archive")
    parser.add_option("-P", "--project_base_dir", dest="project_base_dir",
                      default="/bubo/proj/a2010002/projects")

    parser.add_option("-a", "--flowcell_alias", dest="fc_alias")
    parser.add_option("-y", "--project_desc", dest="project_desc")
    parser.add_option("-l", "--lanes", dest="lanes")

    parser.add_option("-i", "--only_install_fastq", dest="only_fastq", action="store_true",
                      default=False)
    parser.add_option("-f", "--only_install_run_info", dest="only_run_info", action="store_true",
                      default=False)
    parser.add_option("-r", "--run_info_name", dest="run_info_name",
                      default="run_info.yaml")
    # parser.add_option("-r", "--make_delivery_report", dest="report", action="store_true",
    #                   default=False)
    parser.add_option("-c", "--copy", dest="copy", action="store_true",
                      default=False)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    kwargs = dict(
        fc_alias = options.fc_alias,
        project_desc = options.project_desc,
        lanes = options.lanes
        )
    main(*args, **kwargs)
