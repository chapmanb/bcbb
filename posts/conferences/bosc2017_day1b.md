I'm at the
[2017 Bioinformatics Open Source Conference (BOSC)](http://www.open-bio.org/wiki/BOSC_2017)
in Prague. These are my notes from the day 1 afternoon sessions on developer
tools, reproducibility, data science and visualization. It also includes notes
from the keynote by Madeleine Ball on open human genome data.

BOSC is a two day annual conference organized by the
[Open Bioinformatics Foundation](https://www.open-bio.org/wiki/Main_Page). It's
an open community interested in building reusable and interoperable tools for
doing biological research.
All video from these presentations will be available from the
[schedule page](http://www.open-bio.org/wiki/BOSC_2017_Schedule)
and slides from the [f1000 channel](https://f1000research.com/collections/bosc).

# Developer tools / reproducibility

## MultiQC: Visualising results from common bioinformatics tools
_Phil Ewels_

Phil works at SciLifeLab and is part of National Genomics Infrastructure in
Sweden. Movitated by moving from command line based QC done on their many
samples. Developed [MultiQC](http://multiqc.info/), which collects multiple
samples and multiple QC methods into a single report. It provides an at a glance
overview of many metrics and is a brilliant way to present QC. bcbio integrates
MultiQC and it's been a huge improvement. It can scale to thousands of samples
by having static plots instead of interactive; large tables turn into plots.
On the backend MultiQC has great configuration and it is easy to plug in new
tools. At the Codefest, lots of new users got involved and contributed
significant improvements to MultiQC in two days. Wow, there is a plugin for
MultiQC to interact with ClarifyLIMS: pulls in sample metadata and project
information from the LIMS using the API. Coming soon, a new project called
MegaQC that collects results from many MultiQC runs and visualize trends, uses a
local server with database and website.

## NGL – a molecular graphics library for the web
_Alexander S Rose_

[NGL](https://github.com/arose/ngl) is a WebGL based viewer for molecular
graphics. Easy to inject directly into web pages, can use javascript to tweak
and do calculations on the client side. Handles a wide variety of molecular
formats and tons of representations.

## GRAPHSPACE: Stimulating interdisciplinary collaborations in network biology
_Aditya Bharadwaj_

[GRAPHSPACE](http://graphspace.org) is a collaborative way to share network results across multiple
groups. Motivation: due to the complexity it is currently quite hard to do.
GRAPHSPACE allows sharing with collaborators or openly, 300+ users sharing
already so really filling a niche in the biological community. Nice integration
with Cytoscape web. Really solid project wishing I knew more about networks.


## Efficient detection of well-hopping duplicate reads on Illumina patterned
flowcells
_Tim Booth_

Describes how Illumina patterned flowcells work, in terms of technical problems
that can happen during sequencing. Well hopping involves sequences leaking into
surrounding cells. This creates technical duplicates which over-represents the
amount of sequence for RNA-seq or variant calling. Standard approach is Picard
mark duplicates after alignment. Their solution is a quicker and works directly
from bcl files. [Code available on GitHub](http://github.com/EdinburghGenomics/well_duplicates)

## An ensemble approach for gene set testing analysis with reporting capabilities
_Monther Alhamdoosh_

Gene set testing: run an experiment, get a bunch of outputs, then try to
organize into groups that provide biological meaning. Different methods end up
with different gene sets -- hard to identify what is correct, danger of trying
multiple methods until you get something that "makes sense." Ran 12 different
methods and then compared them.
[EGSEA (pronounce eg-zee)](http://www.bioconductor.org/packages/release/bioc/html/EGSEA.html)
supports 12 methods and runs some limma magic to combine into a final gene set
combining the information from all of them. It's a R package available from
bioconductor.

## OpenMS 2.0: a flexible open-source software platform for mass spectrometry data analysis
_Timo Sachsenberg_

Main high throughput method in proteomics/metabolomics in LC-MS, which separates
and quantifies analytes. [OpenMS](http://www.openms.de/) is an open source
framework to compute on mass spec output. Everything open source and well
tested, been developed for 10+ years. Great example of long running successful
project. Has a workflow, tool and library layer to access -- run something
pre-built, combine tools however you want, or code against the API directly. 185
tools that cover almost everything you'd want to do with mass spec.

## Interoperable, collaborative multi-platform variant calling with bcbio
_Brad Chapman_

I talked about how the Common Workflow Language enables [bcbio](https://bcb.io)
to run on multiple platforms.

## Gene Set Variation Analysis in cBioPortal
_Kees van Bochove_

[cBioPortal](http://www.cbioportal.org/) is a cancer genomics tool, currently
supported and developed by [The Hyve](http://thehyve.nl/). GSVA is a method to
turn genes into gene sets for analysis. This is fully integrated into
cBioPortal, supporting queries and heatmaps. cBioPortal has full sets of open
public datasets like TCGA data.

# Data Science & Visualization

## The backbone of research reproducibility: sustainable and flexible tool deployment
_Björn Grüning_

How can we make a FAIR (Find, Accessible, Interoperable, Reusable) set of tool
dependencies. Workflow -- development of code, packaging, then deploying. There
was an explosion of package managers: language specific (R, Python, Perl..) and
system specific (deb, RPM) plus higher level containers (Docker). The current
biology community solution to avoid this:
[conda](https://conda.io/docs/intro.html) a cross platform, language agnostic
installation method. Demonstrates how to build a tool with it.
[Bioconda](https://bioconda.github.io/) is a community building scientific
tools. It has an automated, unified build environment and integrates with other
conda communities like [conda-forge](https://conda-forge.org/). Huge incredible
community maintaining tools. Also autogenerates containers from conda recipes,
making Docker and Singularity containers available. Everything happens
automatically via TravisCI. Works with
[Biocontainers](https://biocontainers.pro/) which provides standards for running
Docker tools.

## Reproducible and user-controlled software management in HPC with GNU Guix
_Ricardo Wurmus_ and _Pjotr Prins_

Ricardo talks about the difference between sysadmins and users. Users have
ad-hoc, volatile and provides no backups/updates. If you need 100%
reproducibility it's hard to isolate from the system even with bioconda,
easybuild. They're not properly abstracted. Docker/containers provide some
solution but hard to compose. [GNU Guix](https://www.gnu.org/software/guix/) is
a functional packaging manager which always provides 100% reproducibility.
Level of abstraction matters, reproducible and safe experimentation, Guix makes
sharing easy, ways to use Guix without root. Guix installs conda.

## A Ubiquitous Approach to Reproducible Bioinformatics across Computational Platforms
_John Chilton_

John follows up on the bioconda discussion with how it enables computational
platforms. HPC-ready, cross-platform, easy to manage/maintain multiple versions
of the same tool. [Galaxy](http://usegalaxy.org) provides the ability to have
explicit containers using Biocontainers. Enables cool things like running
Kubernetes and Singularity. cwltool now has configurable dependency resolution.
Shows example of an amazing workflow porting: an 11 tool ChIP-seq workflow that
got ported from AWS + Docker to HPC without Docker.

## Revitalizing a classic bioinformatics tool using modern technologies: the case of the Cytoscape Project
_Keiichiro Ono_

[Cytoscape](http://www.cytoscape.org/) is a 15 year old project, implemented in
Java. Today, v3.5 released and community sill growing, 300+ associated apps.
It's not an ecosystem of tools. But web is the new platform for platform, how
can it move there without getting rid of the existing ecosystem. Now has a REST
API and CyComponents, which are react-based reusable UI/dataviz components.
Pilot project building off this: CyIndex 2.0 -- a repository for biological
networks.

## The SPOT ontology toolkit : semantics as a service
_Olga Vrousgou_

[SPOT is the Sample, Phenotype and Ontologies team](http://www.ebi.ac.uk/about/spot-team)
at EBI trying to help you annotate your data with ontologies. Ontologies help
with data integration. ZOOMA provides automated mapping to onologies. When it's
wrong, OxO cross reference ontologies. The OLS (ontology lookup service) lets
you find terms when cross referencing doesn't help. Webulous lets you create
ontology terms when non exist. All [tools on GitHub](https://github.com/EBISPOT).

## Biopython Project Update 2017
_Christian Brueffer_

[Biopython](http://biopython.org/) project is an original Open Bioinformatics
foundation project from 1999. Walk through all of the new releases and work this
year. Moved to the BSD license, which is helpful since the Biopython license was
custom and essentially BSD. Three new releases this year with lots of new
contributors. Working on a modernized build process with python wheel format.
Lots of use of continuous integration: TravisCI, Appveyor for windows testing,

# Keynote

## Open Sourcing Ourselves
_Madeleine Ball_

Madeleine starts by describing her background: fixed DNA representation on
Wikipedia, one laptop per child project, more Wikipedia open work sharing
knowledge. At the same time, working in George Church's lab. She got involved
with the [Personal Genome Project](http://www.personalgenomes.org/) -- open
sequence data plus open and extensive phenotype data. Director of Research for
several years, worked a lot on ethical and legal issues. Had full engagements
with participants including returning personal genome. Most genomes we share are
not available to participants. Why don't we do this as scientists? PGP had
amazing communication with participants: community forum, conference plus
discussion directly with Madeleine. Contrast with traditional research projects.

Scary date: how worried should you be about sharing your genome data? Some
possibilities: health issues, ancestry, re-identification with Y chromosomes.
Not a lot of easy answers but lessons to learn. Genomes taught more general
lessons. Uses example of location data of her life in NYC -- this also
re-identifies you and could be potentially problematic. Not genome specific.

PGP consent -- very difficult to do correctly. Had an extensive consent process
with a quiz, which makes people actually understand what they are consenting
too. In the end: people like to share and like to be asked. They have >3,500
individual with public genome and genotyping data. There is always data that
goes too far for you: examples of location, Google Search data.

The current research model: original data holders are universal gatekeepers.
Can't talk to the original individuals either without okay/consent/payment to
the original data holder.

Other approaches: what if everyone was in the same place and could hear you?
What is it was easy to return data? What is users could contribute back? Good
example of adding a HealthKit app for OpenHumans. Idea: have all the data
available with individuals and allow them to share, make it available as they
want. This allows you share data across studies over time.

Time for a concrete example. Dana Lewis and [OpenAPS](https://openaps.org/)
project to monitor her Type I diabetic Artifical Pancreas System. NightScout is
software that monitors nighttime blood sugar to prevent it from being too low
because of 1:20 chance of dying in sleep as Type I diabetic.

OpenHumans offers $5000 grants to do work with Open Humans data and also has a
job available for a programmer. It's a great chance to contribute to open,
innovative work.
