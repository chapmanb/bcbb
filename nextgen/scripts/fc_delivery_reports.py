#!/usr/bin/env python
"""Make delivery notes for a flowcell

Usage:
     fc_delivery_reports.py <flowcell id> 
                       [--archive_dir=<archive directory>
                        --analysis_dir=<analysis directory>]

Given a flowcell id, make delivery reports for all projects on that flowcell.
This script relies on Illumina data being delivered to an archive directory, 
and that the bcbb pipeline has been run in the analysis directory.

The script loads the run_info.yaml and generates a report for each project.

Options:

  -a, --archive_dir=<archive dir>          The illumina data archive root in which 
                                           flowcells are found
  -b, --analysis_dir=<analysis directory>  The directory where bcbb analyses are found
                                           in FLOWCELL_ID directories
  -n, --dry_run                            Don't do anything, just list what will happen
"""

import os
import sys
from optparse import OptionParser

import yaml
import glob
import re
from mako.template import Template
from mako.lookup import TemplateLookup

from bcbio.log import create_log_handler
from bcbio.pipeline import log
import bcbio.templates.mako2rst as m2r
from texttable import Texttable

TEMPLATE="""\
Delivery report for ${project_id}
=================================

Delivery
--------

The clustering was performed on a cBot cluster generation system using
a HiSeq paired-end read cluster generation kit according to the
manufacturer's instructions. The samples were sequenced on an Illumina
HiSeq 2000 as paired-end reads to 100 bp. All lanes were spiked
with 1% phiX control library, except for lane 8, which has 2% phiX.
The sequencing runs were performed according to the manufacturer's
instructions. Base conversion was done using Illuminas OLB v1.9.

Comment
--------

General information
-------------------

${infotable}

${lanetable}

The sequence files are named after the following scheme:
lane_runname_sample_1(2).fastq, where the 1 or 2 represents the first
(forward) and the second (reverse) read in a paired-end run. Single
end runs will have one the first read. The files only contain
sequences that have passed Illuminas Chastity filter.


Run information
---------------

Summary read 1
~~~~~~~~~~~~~~

${read1table}

Summary read 2
~~~~~~~~~~~~~~

${read1table}

QC plots
~~~~~~~~

${qcplots}

${qc30plots}

${errorrate}
"""

def main(flowcell_id, archive_dir, analysis_dir):
    print " ".join([flowcell_id, archive_dir, analysis_dir])
    fp = os.path.join(archive_dir, flowcell_id, "run_info.yaml")
    with open(fp) as in_handle:
        run_info = yaml.load(in_handle)
    project_ids = dict()
    for lane in run_info:
        (l, id) = [x.strip() for x in lane['description'].split(",")]
        if project_ids.has_key(id):
            project_ids[id].append(lane)
        else:
            project_ids[id] = [lane]

    sphinx_defs = []
    for k in project_ids.keys():
        lanes = [x['lane'] for x in project_ids[k]]
        log.info("saw project %s in lanes %s" %( k, ", ".join(lanes)))
        sphinx_defs.append("('%s', '%s_delivery.tex', 'Delivery note', u'Scilife', 'manual'),\n"  % (k, k))
        projectfile = "%s.mako" % (k)
        fp = open(projectfile, "w")
        fp.write(TEMPLATE)
        fp.close()
        mylookup = TemplateLookup(directories=['./'])
        tmpl = Template(filename=projectfile, lookup=mylookup)
        proj_conf = {
            'id' : k,
            'lanes' : project_ids[k],
            'archive_dir' : archive_dir, 
            'analysis_dir' : analysis_dir,
            'flowcell' : flowcell_id,
            }
        d = generate_report(proj_conf)
        rstfile = "%s.rst" % (k)
        fp = open(rstfile, "w")
        fp.write(tmpl.render(**d))
        fp.close()

    sphinxconf = os.path.join(os.getcwd(), "conf.py")
    if not os.path.exists(sphinxconf):
        log.warn("no sphinx configuration file conf.py found: you have to edit conf.py yourself!")
    else:
        fp = open(sphinxconf)
        lines = fp.readlines()
        fp.close()
        sdout = []
        modify_conf = False
        for sd in sphinx_defs:
            if not sd in lines:
                sdout.append(sd)
                modify_conf = True
        if modify_conf:
            i = lines.index("latex_documents = [\n")
            newconf = lines[:i+3] + sdout + lines[i+3:]
            fp = open("conf.py", "w")
            fp.write("".join(newconf))
            fp.close()


def generate_report(proj_conf):
    d = { 
        'project_id' : proj_conf['id'],
        'infotable' : "",
        'lanetable' : "",
        'read1table': "",
        'read2table': "",
        'qcplots': "",
        'qc30plots': "",
        'errorrate': "",
        }

    ## General info table
    tab = Texttable()
    tab.add_row(["Project id", proj_conf['id']])
    tab.add_rows([["Run name:", proj_conf['flowcell']],
                  ["Uppnex project", ""]])
    d.update(infotable=tab.draw())
    
    ## Lane table
    tab = Texttable()
    tab.add_row(["Lane", "Sample(s)", "Conc. (pM)"])
    for l in proj_conf['lanes']:
        samples = []
        for mp in l['multiplex']:
            samples.append(mp['name'])
        tab.add_row([l['lane'], ", ".join(samples), ""])
    d.update(lanetable=tab.draw())
                
    ## qcplots
    byCycleDir = os.path.join(proj_conf['archive_dir'], proj_conf['flowcell'], "Data", "reports", "ByCycle")
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "QScore_L%s.png" % (l['lane']))), width="100%"))
    d.update(qcplots= "\n".join(res))

    ## qc30plots
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "NumGT30_L%s.png" % (l['lane']))), width="100%"))
    d.update(qc30plots= "\n".join(res))

    ## qcplots
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "ErrRate_L%s.png" % (l['lane']))), width="100%"))
    d.update(errorrate= "\n".join(res))
                
    return d

if __name__ == "__main__":
    usage = """
    fc_delivery_reports.py <flowcell id>
                           [--archive_dir=<archive directory> 
                            --analysis_dir=<analysis directory>]

    For more extensive help type fc_delivery_reports.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--archive_dir", dest="archive_dir", default="/bubo/proj/a2010002/archive")
    parser.add_option("-b", "--analysis_dir", dest="analysis_dir", default="/bubo/proj/a2010002/nobackup/romanvg")
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(
        archive_dir = os.path.normpath(options.archive_dir),
        analysis_dir = os.path.normpath(options.analysis_dir)
        )
    main(*args, **kwargs)
