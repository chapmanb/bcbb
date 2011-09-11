<%!
import os
import sys
import glob
import bcbio.templates.mako2rst as m2r
%>

Project
=======

Contents:

.. toctree::
   :maxdepth: 2

Analysis
=============

The Science for Life Laboratory bcbio pipeline 
has been used to map sequences to a reference
sequence. Variant calling has been performed
bla bla bla

Analysis settings
-----------------

${m2r.program_info(proj_conf)}

Mapping results
==============================



Targeted resequencing QC
================================

% if use_picard_metrics

Picard metrics
--------------

Align metrics
~~~~~~~~~~~~~

<%!
for fc in flowcell_ids:
    align_files = glob.glob(os.path.join(intermediate_dir, fc, "*", "*.align_metrics"))
    metrics = dict()
    for af in align_files:
        f = open(af)
        tmp = f.readlines()
        close(f)
        metrics[fc] = tmp
%>
% for fc in flowcell_ids:
${metrics[fc]}
% endfor


Duplication metrics
~~~~~~~~~~~~~

GC metrics
~~~~~~~~~~~~~

Insert metrics
~~~~~~~~~~~~~

SNP filter metrics
~~~~~~~~~~~~~



% endif


% if use_TEQC

TEQC
-----


% endif
