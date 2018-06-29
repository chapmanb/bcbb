I'm at the
[2018 Bioinformatics Community Conference](https://gccbosc2018.sched.com/) in
Portland, Oregon. This is a collaboration between the [Galaxy Community
Conference](https://usegalaxy.org/) and the [Bioinformatics Open Source
Conference (BOSC)](https://www.open-bio.org/wiki/BOSC). These are my notes from
the second day (June 28, 2018).

I'd encourage you to fill out our
[Community Conference Survey](http://bit.ly/gccbosc2018-feedback), whether or not you
were able to physically attend GCC and BOSC this year. This is the first joint
Bioinformatics Community Conference with GCC and an experiment to help improve
what we can offer the community. As a result of this collaboration with GCC, we
now have training days, a more affordable conference, childcare, and up to four
days of collaboration with lots of space to welcome everyone. We have over 180
people signed up for the [CollaborationFest](https://galaxyproject.org/events/gccbosc2018/collaboration/)
with a lot of [great projects](http://bit.ly/cofest2018-ideas), and I'm excited
to see what everyone builds and improves on over the next few days.

We need feedback from the community about what you think about these changes,
and want to learn if co-locating BOSC with GCC, rather than ISMB, helped or hindered
your ability to attend. Please fill out the
[Survey](http://bit.ly/gccbosc2018-feedback); we'd love to hear your thoughts.

# Developer Tools and Libraries

## iMADS: A sustainable software collaboration for predicting transcription factor binding specificity
*Dan Leehr*

Dan discussing their software engineering approach at [Duke GCB](https://github.com/Duke-GCB/), specifically a
collaboration for predicting transcription factor binding based on microarrays.
iMADS collaboration: researcher with science methods, datasets and some initial
python code. Informatics wanted to focus on sustainability, building a minimal
viable product and continue to iterate on it: learn from users using the product
as soon as possible. Milestone 1: generate whole genome predictions and publish
as browser tracks: [Predict TF
Binding](https://github.com/Duke-GCB/Predict-TF-Binding) command line. Milestone
2: run predictions in parallel splitting by chromosomes using a [CWL
workflow](https://github.com/Duke-GCB/TrackHubGenerator); entirely reproducible
and maintainable. Next step was building an online databases with scores and
predictions: database engineering challenges around storing sequence data. [Web
site for visualization and download](https://imads.genome.duke.edu/) and
[website code](https://github.com/Duke-GCB/iMADS). Provides on-demand
predictions using CWL with dockerized code. Wanted to document the process of
building and using this, providing infrastructural knowledge for doing again.
Automated deployment with Ansible and Docker magic.

## The Funnel Task Execution Server
*Alexander Buchanan*

Alex [talking about](http://tinyurl.com/funnel-slides) [the Funnel execution
server](https://github.com/ohsu-comp-bio/funnel). Motivation: create a task,
submit to analysis sysem, assigns to worker which stages inputs, runs tasks, and
uploads outputs. Have status server querying how tasks are doing. Want to take
same task and move over multiple systems. To solve this problem, GA4GH created
task execution schema (TES) as schema and REST API for running tasks. Funnel is
the reference server for this API and contains a large number of plugins for
different backends. Funnel provides a dashboard for monitoring tasks, collecting
CPU and RAM profiles. Has a nice command line interface for creating, managing
and canceling tasks. Funnel has a lot of different backends for databases, file
storage, compute and event management. Simple task API: create, get, list,
cancel. The task message in JSON has metadata, resource requests,
inputs/outputs, commands to run. After the task is run, get ID, status, logs,
uploaded files, and details from the run. Built an ecosystem around this using
JSON + rest, all generalized to TES (rather than Funnel specifically). Bunny and
Cromwell run TES backends using funnel. Galaxy has an experimental task runner
for Galaxy, planning to improve at CoFest next two days. Experimental workflow
engine: [Gaia](https://github.com/bmeg/gaia).

## Distributed execution of bioinformatics tools on Apache Spark with ADAM and Cannoli,
*Michael Heuer*

Michael works in the AMPlab on [ADAM](github.com/bigdatagenomes/adam),
available via homebrew, conda and Docker. Focus is on using Spark for analysis:
core ADAM APIs sit on top of Spark and then enable applications. Talking about
cannoli application today. ADAM transforms the data into efficient distributed
resource: Avro in Parquet as RDDs, also available as Spark Datasets and DataFrames
efficiently partitioned by genomic region.
[cannoli](https://github.com/bigdatagenomics/cannoli) provides a pipe API which
allows streaming in common formats: example of streaming bwa with nice interface.

## Extended Extraction Transform Load: A novel framework for batch jobs on cloud computing resources
*Hiromu Ochiai*

Hiromu talking about genome analysis approaches on cloud resources. Approaches
to using the cloud: building a cluster on cloud (Galaxy, cnf-cluster) -- helpful
since we're used to that approach, but requires persistent static resources so
hard to align with making use of native cloud approaches. Alternative approach:
on demand ETL (Extract, Transform, Load). Tool called [awsub](https://github.com/otiai10/awsub), which starts
instances and pulls from s3 buckets, along with common reference data, executes
and cleans up. Problems with ETL is that the common reference data is large and
requires lots of unnecessary transfer. Approach, use Extended ETL (ExTL) to
pre-fetch and share reference data.

## Plugging Docker-based visualizations into Django with django_docker_engine
*Chuck McCallum*

Working on [Refinery Platform](http://www.refinery-platform.org/) which wraps
Galaxy workflows and feeds into visualizations. Developmd
[django_docker_engine](https://github.com/refinery-platform/django_docker_engine)
which does the work of routing and sending to containers. Practically the
lifecycle is:, you have containers, inputs, visualizations, initialize, use, and
then remove container. Works well if stateless and tied to multiple visualizations.

## The Arachne Graph Database Server
*Kyle Ellrott*

The [Arachne graph database](https://github.com/bmeg/arachne) deals with
property graphs. These are great ways to store dense and linked data: easier to
create long linkage paths vs SQL. Create a BioMedical Evidence graph with
information from places like CiViC and GO. Arachne provides a graph query engine
with multiple backends: Mongo, Postgres, embedded store. Show nice example of
getting RNA and mutation data from cBioPortal. Nice approach handling the
specific kind of tough data you need to deal with when integrating from multiple
data sources.

## Owlery: An easily deployable web service for making reasoning queries over OWL ontologies web-native
*Hilmar Lapp*

Hilmar presenting work done by [Jim Balhoff](https://github.com/balhoff) on
making use of ontologies for querying.
[Owlery](https://github.com/phenoscape/owlery) allows linking ontologies across
multiple domains via semantic reasoning. Allows inference of data implied by
observations. There are a ton of barriers to using these in web applications:
written in low level languages, not designed to communicate externally. So as a
result often use pre-resaoned data, or else need to use a client/server
approach. OWLery provides web service API with JSON-LD and RDF, making it easy
to deploy and integrate. Provides pre-built Docker container for usage with
configuration and linkage to knowledge bases to use.

## taxa: taxonomic data standards and methods for R and Python
*Scott Chamberlain*

[rOpenSci](https://ropensci.org/) is a non-profit developing tools in R for
improving open science. Goal is to deal with taxonomic names: challenges of
heirarchy, multiple sources, often linked to other data by multiple IDs.
[taxa](https://github.com/ropensci/taxa) provides a flexible package for dealing
with taxonomies in R. Nice examples of allowing queries to subset data and taxon
together with a dplyr like interface. [pytaxa](https://github.com/sckott/pytaxa)
is a python version with similar interface, aiming to allow similar
manipulations and transfers between lanuages.

# Workflows

## CWLProv - Interoperable retrospective provenance capture and its challenges
*Farah Zaib Khan*

[Farah's talk](https://slides.com/farahzkhan/cwlprov) focuses on provenance
within workflows. Shows the CWL maintained list of workflows: 215 current
entries. Divide workflows into 3 categories: domain specific pre-built
pipelines, GUI based workbenches, Standardized declarative approaches.
Provenance: information about entities, activities and people involved in
producing data; used to assess trust and re-use of workflows. Provenance
involves everything to replicate and understand the experiment. Different levels
of provenance: level 0 -- nothing about the environment, just the workflow;
level 1 -- content addressible data, executable workflow; level 2 -- multiple
provenance logs, execution details of steps; level 3 -- domain specific
annotations of data. Idea is to combine standards (PROV-model, wfprov and wfdesc
ontology) with CWL into a [ResearchObject](http://www.researchobject.org/). All
standards are interoperable, open source, community driven and domain neutral.
The research object made up of data, snapshots of tool specification files,
workflow with input object and data paths, metadata provenance about the
workflow run. Reference implementation using cwltool. CWLprov helps you
effectively share and re-use your workflow.

## Running portable workflow and container specifications at production scale in the cloud: strategy & best practices, 
*Geet Duggal*

Geet works as DNAnexus but providing general talk about what they've learned
that applies to any approach running on the cloud. Seeing lots of interest in
using CWL and WDL. How can we best write these workflows to allow them to run at
production scale? Scale: tens of thousands of genomes, hundreds of thousands of
exomes. Representing workflows in mostly solved problem with CWL/WDL + Docker.
But what is the experience for executing them when you move between local mode,
cluster HPC system and cloud? Different expectations and limitations depending
on system. Three best practices for writing workflows make it easy to work
across these systems: 1. group together related tasks into computational unit,
but maximize reusability/granularity as much as possible, so tradeoff. 2.
Provide configurability and parameters with the grouped tasks; think less about
stages and more at the high level; differentiate workflow from its steps, think
of utility on its own. 3. Take advantage of CWL and Docker features: filesystem
encapsulation, task/tool definitions, dependencies and conditional execution:
treat workflow-level processing as tasks. Desirable features of an execution
environment: full provenance tracking using Git-like file hashes, automatic
restarts and reusability, versioning publishing and collaboration. DREAM
challenges have helped work out best ways to run CWL/WDL on DNAnexus.

## The GA4GH/DREAM Workflow Execution Challenge
*James Eddy*

James talking about [GA4GH Workflow Execution
Challenge](https://www.synapse.org/#!Synapse:syn8507133/wiki/415976).
Movitation: making tools interoperable which requires standards and APIs plus
implementations. Important question: do these implementations work and give us
the reproducibility we need. DREAM challenges are community competitions; this
one focused on collecting workflows and ensuring they run. Needed to build new
infrastructure in [Synapse](https://www.synapse.org) to manage workflows, data
and checkers. Documented methods, workflows + checking of expected results.

## Bespin: An open source system to run reproducible computational workflows on cloud infrastructure
*John Bradley*

[Bespin](https://github.com/Duke-GCB/bespin) allows domain researchers to
run workflows: CWL bioinformatics workflows + services to run workflows. Provide
generalized methods and references from data stored in CWL representation. Can
automatically generate these. Provides a website, jobs API manages the runs with
a job control backbone to actually allocate resources + data. Runs workflow and
includes a reproducibility kit, outputs, and methods. Initial idea was to
provide full automation, but decided to add review by a bioinformatics expert to
sanity check results and approaches.

## CWL-Airflow pipeline manager as a backend for BioWardrobe data analysis platform
*Andrey Kartashov*

[BioWardrobe](https://github.com/Barski-lab/biowardrobe) is a platform for
analyzing biological data. Enabled scientists to run analyses but ran into
issues with scaling, restarting, versioning and recording metadata. Moved to
using CWL. Looking for pipeline manager to run and utilize, picked Airflow.
[CWL Airflow](https://github.com/Barski-lab/cwl-airflow) integrates CWL with
Airflow. Nice display of dependency graphs, tracking of runs.

## Scaling bioinformatics analysis using Nextflow and AWS
*Francesco Strozzi*

Describing experience adopting Nextflow at [Enterome](http://www.enterome.fr/).
Nextflow has lots of executors: Slurm, Kubernetes, AWS (GCP coming). Uses AWS
Batch with Docker. Automates creation of cluster and handles spot instances.
Computing cost improvements with this and per minute billing. Required minimal
setup on AWS and nice scaling between laptop development and deployment on AWS.
Uses S3 for filesharing which does not scale well in the current implementation.

## One Week to 1,000 Whole Genomes with Open Source: Arvados, CWL, and bcbio
*Peter Amstutz*

Peter describing using [Arvados](https://arvados.org/), developed by
[Veritas](https://www.veritasgenetics.com/) to analyze 1000 whole genomes for
use as controls in a study. Infrastructure challenges. 1. Investing 70Tb of data;
100-way parallel data transfer using elastic compute nodes. Finished in two
hours on AWS instead of 8 days. 2. Organizing the data: use a CWL workflow to
orchestrate data transfer and automate checksums which identified human errors
and swaps. 3. Used bcbio to run the analysis with Sentieon HaplotypeCaller
implementation. At the peak using 8000 cores. 4. Delivering the data: use
parallel transfer with datestamps. 5. Clean up and removal of intermediate data
so it gets cleaned. Records of run kept so have provenance but don't need to
keep files around.

# Galaxy

## Galaxy Community Update
*Dan Blankenberg, Jeremy Goecks, Anton Nekrutenko, James Taylor*

Presenting an overview of the current state of Galaxy, reference [recent summary
paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky379/5001157).
Emphasizes great advancements from the community. Galaxy emphasis on Genomics at
Scale: research is diverse -- example of huge number of individual research
grants from NIH. Makes for lots of dimensions of scale: datasize size, number of
datasets, data/tool heterogeneity. Want to make Galaxy universally usable,
working for folks with smaller communities but important questions. Want to
handle great data where there is not great accessibility. What makes Galaxy
special? Philosophy, accessibility, extensibility, platform mindset and user
engagement. Working on worldwide usage: Galaxy main, Galaxy Europe, Galaxy
Australia. Galaxy training: centered in Freiberg but spread over the world.
Priorities: refocusing on the end user experience. Galaxy UI is unique
advantage: to make it better -- dataset collections, size + heterogeneity,
single page application, multiple connected windows like JupyterHub, drag and
drop, different levels of details with multiple UIs. Planning to take multiple
Galaxy instances and think about federating them; right now each stands
completely alone. Could they interoperate? Authentication/Identification first
problem to solve. Trying to move to Galaxy as a Service. Need to be able to
actual move compute to data where data is all over.

# Late breaking lightning talks

## Common Workflow Language update
*Michael Crusoe*

CWL is a standard for defining workflows, born at CollaborationFest 4 years ago
in Boston. CWL is now a member project of the Software Freedom Conservancy, over
5000 CWL descriptions on GitHub, new CWL user guide. New CWL open source
implementations: Cromwell, REANA (Kubernetes from CERN), CWLEXEC from IBM on LSF
cluster scheduler.

## Code is Science Manifesto
*Yo Yehudi*

Motivation: some scientists don't automatically think scientific code is always
open. We need to be able to see and fix code used in research. Science requires
peer review, science is often computing, software contains errors. Bad code ==
bad science. Session at Mozilla Festival: crowdsourced why people thought code
was closed and what we could do. Idea: produce a manifesto to change culture
around scientific computing. [Code is
Science](https://twitter.com/codeisscience) is the first draft:
[Manifesto](https://codeisscience.github.io/manifesto/): open over closed, code
for future, incorrect code = incorrect science; availability over perfection,
code deserves credit.

## FAIRsharing -- mapping the landscape of databases, standards and data policies
*Massimiliano Izzo*

Mapping a complex and evolving landscape of standards and data. Provides a
web-based and curated portal: collections and recommendations. FAIRsharing is
for researchers and curators who want to be able to find data. Have a
collection/recommendation widget that can be embedded in different websites.

## Scalable computing in Bioconductor: from cores to clusters
*Nitesh Turaga**

Starts by describing map/reduce approach: lapply in R. BiocParallel uses the
same interface: bplapply, BPPARAM parameter determines backend to use for
computation. Example of doing pi approximation, changing between serial and
multicore parallelization. Multiple clusters and backends supported, tons of
customization for specific clusters. Example of doing Salmon pseudoalignment
parallelized on SGE and split by samples. This is the next generation of
parallel approaches in BiocParallel, provides some additional helpers for coding
with it on vectors.

## BioThings Studio
*Sebastian Lelong*

BioThings APIs: mygene.info -- 15 different data sources integrated into merged
JSON database, indexed with elastisearch. Has a REST API for query. Common
backend in Python: [BioThings
Studio](http://github.com/giothings/biothings_studio). To use, create plugin and
parser in GitHub, then can register and generates integration code as a new data
source. Trigger, download and ingest data locally. Studio generates an API and
provides endpoints to query after indexing. Full tutorial:
http://bit.ly/biothings_studio

## Bioinformatics in the Age of AI
*Mike Duncan*

Two talks in one: [MOSES](https://github.com/opencog/moses) a machine learning
tool for addressing issues in genomic data reduction with a genetic algorithm.
Looks at multiple models and rates fitness, evolving answers and avoiding local optima.
[SingluarityNet](https://singularitynet.io/) -- blockchain based network/market
for AI services. The MOSES framework is an example of a potential application on
this network.

# Keynote

## Confound it! Reproducible biology from "omics" data analysis
*Lucia Peixoto*

Peixoto lab uses genomics to understand Autism Spectrum Disorders, specifically
co-occurrence with intellectual disability and sleep impairments. Generate
functional genomic data in mouse, analysis with human WGS data. Everything done
with reproducible methods because they care about doing good science.
Experimental design and challenge of high dimensional data. Good experimental
design minimizes confounders we don't care about and includes controls to
estimate the effect of treatment versus confounders. Case studies: effects of
learning (RNA-seq), effects of sleep deprivation (microarray).

Transcriptomics: most commonly used analysis pipelines for RNA-seq don't remove
confounders other than library size. Recommend RUV as a method to remove
effects. Use normalization to try and remove confounder: differences between
studies, platforms, libraries, day, technician... Poor normalization = more
false positives and negatives. PCA plots reveal confounders. RLE pots reveal
confounders when mean and variance are not similar. Histogram of p-values: if
flat likely have confounders. In neuroscience confounders often dominate the
signal of interest.

What can we do about confounders? Case study on RNA-seq study from mouse for
learning in gene expression. Two learning tasks: fear conditioning, shock mice
when in new environment and can quantify how much they remember this. Leaning
task 2: object location memory, explore area and then move objects on new time.
Get 5 or 6 replicates from two learning methods + controls. Look at standard
pipelines for normalization: UQ, TMM, FPKM, RUV-seq. FKPM and friends do not
change PCA plots, RUV-seq provides real separation. I wish I knew how RUV-seq
compares to DESeq, salmon/kallisto and other more common methods now rather than
length normalization. With pathway analysis, RUV removes encrichment for
incorrect pathways. Second case study, sleep deprivation and recovered sleep.
RUV removing effect of labs and different arrays and recovers most controls.

Epigeomics analysis. Main difference is that you don't have a pre-existing
expecation of where the data will be. Peak calling show low reproducibility
across replicates: well known but unpublished. Need replicates but ENCODE
standard is 2 so a lot of pushback on these. Peak location also confounded.
Brutal reproducibility of H3K9ac across 8 replicates: after 3 or 4 get better
but 1-3 is terrible. Really nice demonstration of importance of replicates. Case
study, how does learning effect chromatic accessibility. RUV normalization also
necessary for differential epigenomic studies. After normalization, learning
regulated regions active during development and associated with splicing.
Compared against ENCODE data expected distributions. DEScan inegrates peak
calling and need at least 3 replicates to use. Still need RUV normalization.
With these two, went from no discovery to providing support. [DEScan2 in latest
Bioconductor](https://bioconductor.org/packages/release/bioc/html/DEScan2.html).

Final considerations: confounders still present even with good experimental
design. Need good design and analyses to avoid issues.
