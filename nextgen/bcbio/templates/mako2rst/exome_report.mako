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
...

Analysis settings
-----------------

${m2r.program_info(proj_conf)}

Mapping results
==============================



Targeted resequencing QC
================================

Picard metrics
--------------

Align metrics
~~~~~~~~~~~~~

${m2r.align_metrics(align_metrics)}

Duplication metrics
~~~~~~~~~~~~~

${m2r.align_metrics(duplication_metrics)}

GC metrics
~~~~~~~~~~~~~

Insert metrics
~~~~~~~~~~~~~

SNP filter metrics
~~~~~~~~~~~~~



TEQC
-----

Summary statistics
~~~~~~~~~~~~~~~~~~~

${m2r.teqc_json(teqc_json)}

Chromosal distribution
~~~~~~~~~~~~~~~~~~

These plots show the distribution of reads over chromosomes.

${m2r.teqc_graphics(teqc_grf)}

Coverage histogram
~~~~~~~~~~~~~~~~~~

The coverage histograms show 

${m2r.teqc_graphics(teqc_grf, which="coverage-hist")}


Coverage uniformity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

${m2r.teqc_graphics(teqc_grf, which="coverage-uniformity")}

Duplicates barplot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

${m2r.teqc_graphics(teqc_grf, which="duplicates-barplot")}

Insert size histogram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

${m2r.teqc_graphics(teqc_grf, which="insert-size-hist")}


Coverage target length
~~~~~~~~~~~~~~~~~~~~~~~

${m2r.teqc_graphics(teqc_grf, which="coverage-targetlength-plot-avgCoverage")}

Coverage target length nReads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

${m2r.teqc_graphics(teqc_grf, which="coverage-targetlength-plot-nReads")}


