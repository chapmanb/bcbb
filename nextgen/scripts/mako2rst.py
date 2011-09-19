#!/usr/bin/env python
"""Convert mako file to rst file for use with Sphinx

Usage:
    mako2rst.py <YAML project config> <mako template file>

"""

import os
import sys
from optparse import OptionParser

import yaml
import glob
import json
from mako.template import Template
from mako.lookup import TemplateLookup

from bcbio.log import create_log_handler
from bcbio.pipeline import log

def main(project_info_yaml, mako_template):
    with open(project_info_yaml) as in_handle:
        project_info = yaml.load(in_handle)
    ## Set the lookup
    _check_conf(project_info)
    mylookup = TemplateLookup(directories=['./', os.path.abspath(mako_template)])
    tmpl = Template(filename=mako_template, lookup=mylookup)
    run_info = get_run_info(project_info)
    d = dict(proj_conf=project_info)

    ## Get TEQC info
    tdata = teqc_graphics(run_info, project_info)
    d.update(teqc_grf = tdata)
    jdata = teqc_json_data(run_info, project_info)
    d.update(teqc_json = jdata)
    ## Get align info
    mdata = picard_metrics(run_info, project_info, "align_metrics")
    d.update(align_metrics=mdata)
    mdata = picard_metrics(run_info, project_info, "duplication_metrics")
    d.update(duplication_metrics=mdata)

    print tmpl.render(**d)


def _check_conf(proj_conf):
    if not proj_conf.has_key('third_party'):
        print >> sys.stderr, "no third_party key: please include relevant values in the third_party section"
    tp = proj_conf['third_party']
    if not tp.has_key('top_dir') or not tp.has_key('intermediate_dir'):
        print >> sys.stderr, "no top_dir or intermediate_dir key defined in YAML config file: need to set project root"
        sys.exit()

def get_run_info(proj_conf):
    flowcells = glob.glob(os.path.join(proj_conf['third_party']['top_dir'], "data", "*"))
    d = dict()
    for fc in flowcells:
        yaml_file = os.path.join(fc, "project_run_info.yaml")
        if not os.path.exists(yaml_file):
            print >> sys.stderr, "no project_run_info.yaml file in %s" % yaml_file
            sys.exit()
        with open(yaml_file) as in_handle:
            run_info = yaml.load(in_handle)
        d[os.path.basename(fc)] = run_info
    return d

def teqc_json_data(run_info, proj_conf):
    jdata = {}
    for fc in run_info.keys():
        jdata[os.path.basename(fc)] = {}
        samples = run_info[fc]
        indir = glob.glob(os.path.join(proj_conf['third_party']['intermediate_dir'], "*" + fc))[0]
        infiles = glob.glob(os.path.join(indir, "*.json"))
        for f in infiles:
            fp = open(f)
            jd = json.load(fp)
            fp.close()
            jdata[fc][os.path.basename(f)] = jd
    return jdata
            

def teqc_graphics(run_info, proj_conf):
    tdata = {}
    for fc in run_info.keys():
        samples = run_info[fc]
        indir = glob.glob(os.path.join(proj_conf['third_party']['intermediate_dir'], "*" + fc))[0]
        png = {'chrom-barplot': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-chrom-barplot.png"))],
               'coverage-hist': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-hist.png"))],
               'coverage-targetlength-plot-avgCoverage': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-targetlength-plot-avgCoverage.png"))],
               'coverage-targetlength-plot-nReads': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-targetlength-plot-nReads.png"))],
               'coverage-uniformity': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-uniformity.png"))],
               'duplicates-barplot': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-duplicates-barplot.png"))],
               'insert-size-hist': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-insert-size-hist.png"))]
                   }
        tdata[fc] = png
    return tdata

def picard_metrics(run_info, proj_conf, ext, match=""):
    """Read picard metrics files for a given extension and match"""
    mdata = {}
    for fc in run_info.keys():
        mfiles = glob.glob(os.path.join(proj_conf['third_party']['intermediate_dir'], fc, "*%s*.%s" % (match, ext)))
        mfdata = {}
        for mf in mfiles:
            f = open(mf)
            tmp = _filter_lines(f.readlines())
            f.close()
            mfdata[os.path.basename(mf)] = tmp
        mdata[fc] = mfdata
    return mdata

def _filter_lines(lines, comment_char="#"):
    lret = []
    for l in lines:
        if not l.startswith(comment_char):
            data = l.split()
            if len(data) > 0:
                lret.append(data)
    return lret

if __name__ == "__main__":
    usage = """
    mako2rst.py <YAML project config> <mako template file>
"""

    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    kwargs = dict(
        )
    main(*args, **kwargs)
