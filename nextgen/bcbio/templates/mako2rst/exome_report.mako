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

Picard metrics
--------------

Align metrics
~~~~~~~~~~~~~

${m2r.align_metrics(align_metrics)}

Duplication metrics
~~~~~~~~~~~~~

GC metrics
~~~~~~~~~~~~~

Insert metrics
~~~~~~~~~~~~~

SNP filter metrics
~~~~~~~~~~~~~



TEQC
-----

Chromosal distribution
~~~~~~~~~~~~~~~~~~

These plots show the distribution of reads over chromosomes.

${m2r.teqc_graphics(teqc_grf)}

Coverage histogram
~~~~~~~~~~~~~~~~~~

Coverage target length
~~~~~~~~~~~~~~~~~~~~~~~


Coverage target length nReads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



