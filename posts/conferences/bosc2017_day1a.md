I'm at the
[2017 Bioinformatics Open Source Conference (BOSC)](http://www.open-bio.org/wiki/BOSC_2017)
in Prague. These are my notes from the day 1 morning session on workflows. It
was an incredible set of talks around the community work done to establish
standards for running workflows, as well as practical implementations that
support them.

BOSC is a two day annual conference organized by the
[Open Bioinformatics Foundation](https://www.open-bio.org/wiki/Main_Page). It's
an open community interested in building reusable and interoperable tools for
doing biological research.
All video from these presentations will be available from the
[schedule page](http://www.open-bio.org/wiki/BOSC_2017_Schedule)
and slides from the [f1000 channel](https://f1000research.com/collections/bosc).

Nomi Harris starts the day off with a welcome to the 18th annual BOSC
conference, describing the goals of the conference. Next year BOSC will partner
with the Galaxy Community Conference (GCC) for
[a bioinformatics community conference](https://gccbosc2018.sched.com/) in
Portland, Oregon June 25-30, 2018.

# BOSC opening talks

## The Open Bioinformatics Foundation
_Hilmar Lapp_

Hilmar describes the goals of the
[Open Bioinformatics Foundation (OBF)](https://www.open-bio.org/wiki/Main_Page).
He polls the audience and 1/3 to 1/2 of attendees are at BOSC for the first
time. We have a great turnout and it's exciting to see so many new people here.
OBF has a
[OBF Travel fellowship program](https://github.com/OBF/obf-docs/blob/master/Travel_fellowships.md)
to help increase diversity, in all aspects, at OBF or other open source
bioinformatics events like BOSC.

## OBF in the Google Summer of Code. Wrapping up 2016 and presenting the 2017
_Kai Blin_

Kai describes OBFs participation in
[Google Summer of Code](https://summerofcode.withgoogle.com/). Google is kind
enough to pay students to work on open source code for the summer. Kai organizes
OBFs efforts to expand supported students beyond OBF's core projects. 7 students
in 2016 on a variety of project and 8 in 2017. This takes an incredible amount
of work from both mentors and students and is a great training process. Already
thinking about 2018 -- check out [the GitHub repository with details](https://github.com/OBF/GSoC).

## Codefest 2017 summary
_Brad Chapman_

I presented the
[work finished](https://github.com/chapmanb/bcbb/blob/master/talks/bosc2017_bcbio_interoperate/chapmanb_bcbio_interoperate.pdf)
at [Codefest 2017](https://www.open-bio.org/wiki/Codefest_2017)

# Workflows

## Rabix Executor: an open-source executor supporting recomputability and interoperability of workflow descriptions 
_Janko Simonovic_

Janko talks about the
[Rabix executor (also known as bunny)](https://github.com/rabix/bunny), an open
source project to run common workflow language. Part of an overall suite of
interoperable parts. They have CWL bindings that integrates the executor with
backends supporting GA4GH TEST, LSF, Slurm. They talk to the queues on the
backend using local, RabbitMQ and ActiveMQ. Using TES, the engine and scheduler
deploy as one component then spawn off executors. The have useful approaches for
handling difficult workflows: subworkflows with look aheads

## Rabix Composer: an open-source integrated development environment for the Common Workflow Language
_Ivan Batic_

[Rabix composer](https://github.com/rabix/composer) is a visual code editor for
the Common Workflow Language. Incredibly feature rich, GUI for creating tools,
inspect the CWL, see your command lines, validation, ability to use test data,
support for editing expressions and a full code editor. There is a graphical
Workflow Editor representation where you can inspect the workflow, inspect
individual steps and inputs/outputs, work through the history and provides
undo/redo as you build the graph, dragging and dropping of files and exports as
SVG. It fully integrates with Seven Bridges, the Cancer Genomics Cloud and
Gavatica (a pediatrics cloud). Next up is work to integrate directly with rabix
executor.

## CWL-svg: an open-source workflow visualization library for the Common Workflow Language
## CWL-ts: an open-source TypeScript library for building developer tools for the Common Workflow Language
_Maja Nedeljkovic_

Maja talks through the libraries used inside Rabix composer that are available
on their own. [CWL-ts (typescript -- editing)](https://github.com/rabix/cwl-ts)
[CWL-svg (visualization)](https://github.com/rabix/cwl-svg). Brilliant
decoupling of concerns and incredibly useful starting points for others working
on building CWL tools

## The GA4GH Tool Registry Service (TRS) and Dockstore - Year One
_Denis Yuen_

Denis talks about an interoperable standard from GA4GH called the
[Tool Registry Service (TRS)](https://github.com/ga4gh/tool-registry-schemas).
Practically, [Dockstore](http://dockstore.org/) is a freely available
implementation you can use now, developed in the
[Pan Cancer Analysis of Whole Genomes project (PCAWG)](https://dcc.icgc.org/pcawg#!%2Fmutations).
Denis ends with a request for people to come to the TRS group

## The GA4GH Task Execution System (TES) and Funnel Server
_Brian O'Connor_

The
[GA4GH task execution service (TES)](https://github.com/ga4gh/task-execution-schemas)
is a way to execute a tool in a variety of environments. It builds off the
registry service (TRS) to provide a way to run a tool in a different
environment. There is a practical implementation called
[Funnel](https://github.com/ohsu-comp-bio/funnel) that you can use right now.
Goal is to build both good standards and tools you can use in production.

## The GA4GH Workflow Execution Service (WES)
_Peter Amstutz_

The
[GA4GH Workflow Execution Service (WES)](https://github.com/ga4gh/workflow-execution-schemas)
builds off of the task execution part to run an entire workflow. So TES is a
step, and WES is the workflow. The idea is to submit a CWL or WDL workflow and a
description of the data and have a way to run it re-usably over multiple
systems. The practical implementation,
[Workflow Service](https://github.com/common-workflow-language/workflow-service)
currently supports Arvados and cwl-runner, work in progress to handle bunny.

## The GA4GH/DREAM Infrastructure Challenges
_Brian D. O'Connor_

Brian summarizes all the standards and tools just discussed and how they're
building blocks to providing implementations. The next goal is to demonstrate
that these tools work on real workflows. The
[GA4GH DREAM infrastructure challenge](https://www.synapse.org/#!Synapse:syn8507133/wiki/415976)
opened yesterday with a bunch of workflows, including bcbio. They need
participants -- take these workflows and run them.

## Workflows interoperability with Nextflow and Common WL
_Kevin Sayers_

[Nextflow](https://www.nextflow.io/) is a domain specific language for
describing bioinformatics workflows and a runner that integrates with multiple
locations. Kevin has worked on [supporting the Common Workflow Language with
Nextflow](). Work in progress that converts tools to nextflow steps and workflows
to nextflow processes. Supports scatter with channels. Working on supporting
expression tools, record inputs/outputs, sub-workflows.

## CWL Viewer: The Common Workflow Language Viewer
_Stian Soiland-Reyes_

Stian describes a problem they identified -- hard to get an overview of a CWL
workflow from a GitHub repo of CWL files to a visual description of a project.
The [Common Workflow Language viewer](https://view.commonwl.org/) will take a
URL and provide a beautiful display including documentation of inputs/outputs
(if folks include them, which they will if made available). Stian also talks
about provenance based [Research Object](http://www.researchobject.org/)
integration with CWL during the Codefest.

## Screw: tools for building reproducible single-cell epigenomics workflows
_Kieran O'Neill_

Keirin folks on perfect reproducibility: need reproducible code and software.
Practically describes an implementation for single cell epigenomic
bisulfite sequencing, focused on cancer research. Screw currently implements one
example data set, plans to work on many public datasets. All available on
GitHub: https://github.com/epigenomics-screw

## BioThings Explorer: Utilizing JSON-LD for Linking Biological APIs to Facilitate Knowledge Discovery
_Jiwen Xin_

Problem identified in translational medicine -- how to find, update and link
data. The key issue is API inteoperability: nice example from Ensembl,
mygene.info, myvariant.info. The Ensembl ID present in 3 different places in the
three schemas with different listing. Makes it a pain to link these together.
JSON-LD resolves this problem by defining key IDs assigning the fields to identifiers.
[BioThings explorer](http://biothings.io/#) uses this to aggregate and organize
an amazing set of data sources. This is an amazing resource.

## Discovery and visualisation of homologous genes and gene families using Galaxy
_Anil S. Thanki_

Gene families: a set of similar genes formed by duplication of a single original
gene. Various tools available to overview gene families but not easy to
investigate structural changes between genes. [Ensembl GeneTrees](http://www.ensembl.org/info/genome/compara/index.html) great tool but
requires a lot of dependencies and hard to setup. Developed a workflow called
[GeneSeqToFamily](https://f1000research.com/posters/5-1571) that implements the
Ensembl pipeline within Galaxy. Provides visualization and evaluation of
outputs. [The workflow is on GitHub](https://github.com/TGAC/earlham-galaxytools/tree/master/workflows/GeneSeqToFamily.)

## YAMP : Yet Another Metagenomic Pipeline
_Alessia Visconti_

[Yet Another Metagenomics pipeline (YAMP)](https://github.com/alesssia/YAMP) is
a nextflow and Docker based workflow to process metagenomic data. First analysis
block is quality
control: de-duplication (clumpify), trimming (BBduk) and decontamination
(BBwrap), FastQC. Second analysis block is community characterization: taxonomic
binning and profiling (metaphlan2), functional annotation.
