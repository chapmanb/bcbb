I'm at GCCBOSC 2018, the
[Bioinformatics Community Conference](https://gccbosc2018.sched.com/) in
Portland, Oregon. This is a collaboration between the [Galaxy Community
Conference](https://usegalaxy.org/) and the [Bioinformatics Open Source
Conference (BOSC)](https://www.open-bio.org/wiki/BOSC).

The goal of the conference is to bring together a community working on
interoperable tools to support biological research. The conference features two
days of training, two conference days of talks and discussion (with some joint
and some parallel sessions), and up to four days of collaboration. It's a
wonderful group of people dedicated to working better together and the first
joint Bioinformatics Community Conference.

These are my notes from the first evening (June 26, 2018) and the first full day
of the meeting (June 27, 2018). (When there were parallel GCC and BOSC sessions,
I attended the BOSC sessions, so my report does not cover the parallel GCC
sessions.)

There are additional[notes on Day
2](https://github.com/chapmanb/bcbb/blob/master/posts/conferences/bosc2018_day2.md).
Also please fill out our [conference survey](http://bit.ly/gccbosc2018-feedback)
so you can have your say in the future of BOSC.

# Keynotes

## Democratizing Data
*Tracy Teal*

[Tracy's
talk](https://gccbosc2018.sched.com/event/EQF7/opening-keynote-tracy-k-teal-democratizing-data)
starts by discussing the positive impact of bioinformatics and current
frustrations. Data and tools are not limiting bioinformatics: it's how and who
does the work. Different framework for software in science, rather than the rest
of our life. Software is not a service; skills and people are integral to the
research process. What can we do to support this? Building local talent:
distributing skills within our community and complement skill sets. Don't try to
hire wonder-woemen, instead share skills and knowledge. Many people have skills
needed, but might just lack the computational skills. Speaks to importance of
training. Challenge: everyone wants training but are busy, so try to train in
the gaps. Needs to be: accessible, approachable (at the right level), aligned
with interests, applicable (use right after they learn it). Practically:
workshops, short courses, MOOCs.

[Software and data carpentry workshops](https://software-carpentry.org/) teach
core skills in 2 day events with trained instructors. 337 workshops in 2017,
7217 learners with global reach (including Antarctica). By democratizing data
skills: answer more questions, more diverse set of researchers, more creativity.
Idea is to create a program that builds talent. General support structure for
services and resources.

Collaboration: how to work together without going in a lot of different
uncoordinated directions. Open code, data and communication. General idea is to
build a community of practice.

Really inspiring talk focusing on how we can support and build each other up to
have a stronger community. Great Q/A discussion to move training forward: how do
we go from introduction to R to machine learning, say. Want to progress the
baseline of what people know so we can move to the next level.

## Sustainable development of scientific open source tools: a view from Jupyter
*Fernando Pérez*

Fernando [focuses on open source lessons from the Jupyter
project](https://gccbosc2018.sched.com/event/EQ3d/keynote-fernando-perez-sustainable-development-of-scientific-open-source-tools-a-view-from-jupyter)..
He started in open source development for four reason: ethical -- openness is
fair; human/social -- openness fosters collaboration; epistemological --
proprietary science is an oxymoron; technical -- fun to work in Python.
Programming is primarily about communication between humans and we're not
scaling our ability to understand code due to fundamental limitations in
understanding and collaboration. IPython started as idea of iterating on code in
console, now a big team project. Jupyter a web based front end mixing code, text
and results with a formalized communication protocol. This now means it's
language agnostic and can get used from any programming language: 100 languages
supported. Lots of UI front ends, including the classic notebook. [JupyterLab](https://github.com/jupyterlab/jupyterlab):
collects all of the lessons into a single tool. Document is a source of
computation, so more than just tools glued together. Nice demo of accessible
large CSV inputs without lag. Nice viewer integration of FASTA pileup display.
Can use data viewers as code, really nice choices for how to substitute
components. Shows markdown editing with live previews, fully computational
environment essentially breaking up monolith in original application.

Incredible amount of work and refactoring in the current version, a really
impressive display of vision and implementation. The other impressive thing is
the number of third party applications developed on the framework that
interoperate seemlessly.

[JupyterHub](https://github.com/jupyterhub/jupyterhub) is a centralized place
for coordinating organizations. Can use shared hardware/compute and code across
an authenticated setup. Nice example of interactive distributed deep learning
run on big supercomputer with live monitoring. Large scale run plus terminal
based interactivity. JupyterHub also provides ability to access to code
somewhere. Then can package this all up in [binder](http://mybinder.org):
Jupyter + GitHub + explicit dependencies, then automatically creates Docker +
Kubernetes magic to make it available from a single click. Reproducible
research: complete software development environment and instructions that
generate the figure.

This all helps make computational hygiene a habit, rather than something you do
at the end of a process. Couple this with training to give people the ability to
get started in this type of reproducible environment from the beginning.
Re-using training materials from Software Carpentry. Teaching a project based
setup that uses standard repos and practices, gives a standard playbook to learn
how to get going with reproducibility.

Challenges in the Jupyter community: better at building communities, inclusion
and sustainability. Industry partners have a responsibility, do not behave in a
extractive way of taking key people from open source projects. This is toxic to
building an inclusive community; not everyone can work for free to demonstrate
value. Lots of hard problems: Jupyter is successful, but sustainability always
in question. Deep academic connections but models at odds with open source
models. These are sources of technology and economic value but no standard
business model. Good read -- [Nadia Eghbal: an alternative ending to the tragedy
of the commons.](https://nadiaeghbal.com/tragedy-of-the-commons).

# Galaxy packaging

## The journey of a team of engineers in learning packaging technology
*Laure Quintric, Patrick Durand, Valentin Marcon*

Laure representing 4 engineeers in 2 different places in France. [Their
talk](https://gccbosc2018.sched.com/event/EYCM/the-journey-of-a-team-of-engineers-in-learning-packaging-technology)
discusses how they provide bioinformatics tools to their scientific communities.
How do they package and manage tools? Example with [FROGS, a metabarcode analysis
pipeline](https://github.com/geraldinepascal/FROGS) in Galaxy. Widely used but
difficult to manage the installation which had a ton of manual steps.
Technologies used to improve: conda/bioconda, Planemo. Successful deployment in
bioconda plus a bunch of lessons learned in terms of development: tests, logging
frameworks.

# Translational and Medical Informatics

## Cloud-based data ingest for the Human Cell Atlas
*Daniel Vaughan*

[Dan
talks](https://gccbosc2018.sched.com/event/Eit5/cloud-based-data-ingest-for-the-human-cell-atlas)
about the [Human Cell Atlas](https://www.humancellatlas.org/), a large scale
collaborative project for mapping human cells. Dan leads the
data ingest and organization at the EBI: painless submissions plus high quality
FAIR data on top of scalable, flexible architecture. For the metadata model, try
to stay flexible: projects, biomaterial, processes, files and protocols. Uses
JSON all the way down so metadata separate from schema and can evolve independently.
Uses spreadsheet conversion for specification and ingest. Uses Kubernetes to
allocate resources. To productionize this, use Terraform/Helm to manage and
moving towards managed Kubernetes platforms. [NeMO Neuroscience Multi-Omic
Archive](https://nemoarchive.org) already running this.

## The Human Cell Atlas Data Coordination Platform
*Brian O'Connor*

[Brian's
talk](https://gccbosc2018.sched.com/event/Eisz/the-human-cell-atlas-data-coordination-platform-an-open-source-cloud-based-system-for-ingesting-storing-analyzing-discovering-and-exploring-single-cell-data)
focuses on the Data Coordation Platform (DCP), an open source, cloud-based
system for ingesting, storing, analyzing, discovering, and exploring single-cell
data, Store sits on top of object
stores with multi-cloud synchronization, unform cloud access, indexing and
search across metadata, data and metadata organized into bundles. Data bundles
have versioning and grouping of files. Data Store has a RESTful API for access
to files and searching. Also has a command line (`pip install hca`). The
analysis component sits on top of DCP with automated analyses. Uses WDL-based
workflows registered in Dockstore with Cromwell as a service on Google Cloud.
Working on a benchmarking environment to compare workflows and compare against
production tools. Encouraging best practices. Discovery component of DCP is a
portal with links to documentation and a data browser to quickly find subsets of
HCA data. Don't want data browser to be the only place to look at data, want to
have multiple platforms. So working on building tertiary analyses provided by
other groups: needs query and handoff of results. There are already preview
datasets available, and APIs already available. [Data
Biosphere](https://medium.com/@benedictpaten/a-data-biosphere-for-biomedical-research-d212bbfae95d)
is a larger idea that there are shared components across multiple projects.
Collaboration: can projects pick up and use components and also use API
standards. All of the code is freely available https:/github.com/HumanCellAtlas
and data.

## Open Humans - connecting, sharing and analyzing personal data that enables community-driven research
*Bastian Greshake Tzovaras*

[Bastian discusses personal data
use](https://gccbosc2018.sched.com/event/Eity/open-humans-connecting-sharing-and-analyzing-personal-data-that-enables-community-driven-research)
in the context of genomics, trying to avoid the manipulation and monetization of
data happening with companies like Facebook. Alternative approach: collect data
directly from individuals, enabling existing communities. A better approach to
engage with diverse communities. Participant-centered research, involved in
designing the study and can more effectively give consent. How can we do this on
the web at scale? [Open Humans](https://www.openhumans.org/) provides a platform
to connect data and people with projects. Projects can be run by anyone, so open
o researchers and citizen scientists. Examples of projects: activity trackers,
easy way for researchers to ask for permissions. Genome exploration: Genevieve
personal genome annotation with Clinvar. Deep connections with Jupyter notebooks
for integration and anslysis, correlating data sources across multiple inputs.
Personal Data exploratory is a location for sharing notebooks.

## Enabling Machine Learning at Scale with “Tiled” Human Genomes,
*Sarah Wait Zaranek*

[Tiling is an approach to break up a genome into shorter
sequences](https://gccbosc2018.sched.com/event/Eiu0/enabling-machine-learning-at-scale-with-tiled-human-genomes).
Advantage is to recast variants into tiles and then represent as an array of
tiles across a set of tiles. This creates a numerical matrix you can use out of
the box with machine learning methods. Lightning is the tool that manages these,
representation + software on Arvados and CWL. Test Study on Personal Genome
Project data: build a SVM to classify eye color and antigen classifiers based on
known data. Proof of concept resolution on tiles related to known genes. Second
test studey on Alzheimer's data, focus on chr19 with ApoE, classifier finding
known risk factor for ApoE4. Currently scaling approaches to full genome.
Takeways: out of the box ML methods with tiled genomes, next steps to extend to
more complex models. Open source tiling code: https://github.com/curoverse/l7g

## JASS: a free and open source software for the joint analysis and interactive analysis of GWAS results
*Hervé Ménager*

[JASS allows analysis and exploration of GWAS
results](https://research.pasteur.fr/fr/event/joint-analysis-of-multiple-phenotypes-using-gwas-summary-statistics-by-hugues-aschard/).
Uses joint analysis to improve statistical power of results. Python package to
combine GWAS summary statistics. Creates Manhattan plots on command line or via
a web interface with zooming and other fancy visualizations.
http://jass.pasteur.fs and also available locally.

## Variant Transforms: Loading VCF files into BigQuery for large scale data analysis
*Asha Rostamianfar*

[GCP Variant
Tranforms](https://github.com/googlegenomics/gcp-variant-transforms) is a tool
for scalably converting VCF files into Google's BigQuery. Once you're in
BigQuery than can do a ton of cool scalable query on top of it.
VariantTransforms uses Apache Beam and handles malformed/incorrect VCFs
including invalid records and headers. Built a scalable VCF validator tool that
lists inconsistencies and problems. Infers a schema with typing based on the
inputs. BigQuery supports standard SQL and shows a cool example of gender
inference across all of the samples. Has native support for parsing annotations
from VEP and snpEff.

# All About Data

## miRTop: An open source community project for the development of a unified format file for miRNA data
*Lorena Pantano Rubino*

Lorena talking about the [miRTop](https://github.com/miRTop/mirtop)
collaboration. Goal is to detect and annotate miRNAs from short read sequencing.
miRNAs involved in gene regulation. miRNAs have a diversity of isomiRs for each
target to a mature mRNA. Tools to manage this, but each developed independently.
Ideally have the steps defining standards for each step in the pipeline and then
multiple implementations for each. MirGeneDB in GFF3 format, defining
attributes. mirtop handles stats, joining GFF files, comparing files and
converting to counts for expression analysis. Using public data sets for
validation with replicates across multiple labs and tools. Counts for
variability of total number of miRNA reads, stratify across different labs for
reproducibility with a single tool. iso_snps very variable due to sequencing
errors. Next step is development of filters to remove non reliable sequences.

## Reproducible big data science: A case study in continuous FAIRness
*Ravi K. Madduri*

Ravi has two components to talk: 1. reproducibility and tools to help make it
happen, 2. Case study of data availability. Echos Fernando's discussion this
morning: very hard to make your work reproducible after the fact. Case study
required a ton of effort. Need to be baking this into our work from the
beginning. Parallel between reproducibility and teaching kids to brush their
teeth. It's easy once it becomes a habit but need good tools and on ramp to
getting into making that possible. Globus Auth: foundational serivce for
authentication and authorization ecosystem. Used across NIH Data Commons for
data sharing and transfer. Two elements of interoperability in the FAIRness
experiment: minimal identifiers for tracking files (Minid) and self-describing
formats for exchange ([BDBag](http://bd2k.ini.usc.edu/tools/bdbag/)). Workspaces
run analyses at scale using Galaxy and Jupyter notebooks. Case study examples
use transcription factor binding predictions across multiple methods. Lessons
learned from proecess: how do we manage data lifecycle and long term storage,
reproducibility requires planning.

## InterMine 2.0: More than fifteen years of open biological data integration
*Yo Yehudi*

[InterMine is a open source data warehouse
system](https://github.com/intermine/intermine) designed to normalize and
integrate multiple data sources. A whole page (34) of different databases
running InterMine for multiple organisms and data types. Technical stack: Java
based front end stored in Postgres backend optimized for reads. 8 language
toolkits to work with data. GUI for visualizations and statistical enrichment.
Now have new [BlueGenes interface](https://github.com/intermine/bluegenes) built
with Clojure and ClojureScript. Community is one of the most important parts of
OpenSource: who and how are you communicating with? Highlight contributors.
Value of Google Summer of Code students in expanding contributors and community.


## GRADitude: A computational tool for the analysis of Grad-seq data
*Silvia Di Giorgio*

Grad-seq separates RNA-RNA and RNA-protein complexes based on molecular weights.
Can see similarities between experiments using clustering.
[GRADitute](https://github.com/konrad/GRADitude) is the tool that provides
plotting and analysis of data. tSNE plots, RNAs and protein combined plots with
correlations. Next steps are to share the data and make it public.

## NIH Data Commons Pilot Phase leverages the cloud to access, analyze, and share FAIR biomedical data
*David Siedzik*

David presenting an overview of the NIH Data Commons. Data Commons is an effort
to address the increase in data that needs to be made available in a FAIR
manner. It's a new way of funding platform development that focuses on
collaboration and interoperability. It's also a development of a framework. 3
pilot datasets: GTEx, TopMed and Alliance of model organism databases. Using a
lot of standards for interoperability and developing community governance and engagement.


## Reproducible data analysis with Snakemake
*Johannes Köster*

Snakemake is a widely used workflow system, available on bioconda. Concise
python-based DSL for defining inputs and outputs. Can write declarative
workflows and expand with full power of Python. Handles local and remote file
inputs. Can use re-usable tool wrappers or CWL tool definitions. Automatically
determines dependencies and creates a DAG that can get parallelized on clusters
and cloud resources. Snakemake also has full conda integration for defining
dependencies and the environment for steps. Also has Singularity integration for
packages and can combine both conda and Singularity for maximum control.

## Improving the Bioinformatics Curriculum
*Jason Williams*

Conspiracy theory on training: bioinformatics curriculum dominated by unseen
forces. Lots of centers of gravity driving training: tech sector, biology domain
knowledge, math/statistics. Focus needs to be on people: brain is the same, no
matter which sector we look at. We accelerate training by teaching computational
skills; 93% of people will be working with large datasets but only 12% rank
themselves as advanced in dealing with this data. Large survey of 1260 people:
behaviors in addressing bioinformatics. 95% believe bioinformatics should be
part of curriculum by only 40% actually do this. 41% of underrepresented
minorities report training barriers. Minority serving institutions different in
terms of training and ability to train: only 23% report integrating
bioinformatics versus 43% at non-minority serving institutions. New faculty
aren't integrating bioinformatics: more training the newer you are, but less
percentage of faculty integrating bioinformatics. How can we address challenges?
Meet faculty where they are: bioinformatics core competencies to teach in
undergrad, incubator to help create plug and play lessons. Opportunities:
community colleges, undergraduate only institutions, target students early in
degrees, underrepresented students: sell skills they can use like data science
that are applicable to bioinformatics. Keys: effectively communicating with people.

