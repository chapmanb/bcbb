I'm at the first NIH Data Commons workshop in Boston at Countway Library of
Medicine, May 30-31 2018. The goal of the NIH Data Commons is to build an
engaged community of connected software and data. [Titus Brown's overview of the
Data Commons](http://ivory.idyll.org/blog/2018-growing-community.html) is a good
place to get started, as is [Rayna Harris' summary of the organizational
strategy of the Data Commons
Consortium](https://medium.com/@raynamharris/open-source-style-community-engagement-for-the-data-commons-pilot-phase-consortium-f959abe7c0c5).

Vivien Bonazzi from the NIH starts off the meeting taking questions from the
audience about what they want to know about the commons. Since it's a new
effort, most of these focus on defining the goals of the commons and what
success means for it. The general trend at NIH is that big projects get funded
and produce data, but then the data and tools don't necessarily interoperate at
the end of the project funding. How can we change this so we build systems that
work from each other and continue over multiple projects. What infrastructure
layers give us the ability to handle data at scale. Hardest path is now to have
flexible open infrastructure; want to do this using the grass roots from the
ground up. Community should be building this, not top down from NIH. This is a
great opportunity to design the systems we want to have with support from NIH.

# Project Lightning talk updates

The first morning, groups presented short lightning talks about ongoing work to
this point.

## Team Copper -- Training and Project Management
/Charles Reid/

Team Copper is responsible for training, reviewing and infrastructure: UC Davis,
Curoverse and Harvard Chan. Planning events and workshops. Generally operating
on GitHub at https://github.com/dcppc (workshops, project management). Open
source centric group: code and project on GitHub, Group Chat on Slack, e-mail on
groups.io and Zoom for teleconference. Goal is to increase communication with
the larger biomedical community.

## Team Argon -- Data organization and transfer via Globus and ecosystem
/Carl Kesselman/

Overall approach: full stack team providing an ecosystem around the Globus
platform. Data control, movement and workflow management. Backend for data
sharing efforts at the NIH. Uses Big Data Bags (BDBags) -- simple mechanisms for
describing data and large collections. Standard ways of aggregating and
transfering data.

## Team Calcium -- Data Biosphere vision
/Cricket Sloan/

Full Stack group with 3 stacks: Broad, UC Santa Cruz and Vanderbilt. Data
Biosphere: modular, community driven, open and community based. Working closely
with GA4GH with APIs for moving around data. UCSC stack: CGP -- Dockstore;
UChicago: Gen3: Windmill, IndexD; Broad: Firecloud, workspace area. Current
work: coordinating between these 3 platforms: workflows from TopMed to
Firecloud; found data from Boardwalk and Windmill to bring to FireCloud and run
analysis there using the Cromwell Workflow Engine.

## Team Carbon -- science use case with data search, access and central APIs
/Paul Avillach/

Proving scientific use case: automatically classify variations associated with
heart diseases across diverse populations. Need data commons to be able to
search, retrieve data, compute on data, and distribute results, all across a
diverse ecosystem. Providing an API for finding, moving and getting data. Have a
simple implementation supporting this. SmartAPI with 10 already registered APIs.

## Team Phosphorus -- multi-cloud data infrastructure middleware using existing APIs
/Claris Castillo/

Full Stack team, establish a science driven group with uses cases. Exploring
different data access methods and APIs. Built a Cross-Cloud Compute:
Mesos/Marathon/Kubernetes/Docker stack. Next step is data management across
multiple clouds using iRODs. Build Middleware services to bring these together:
PIVOT API/Core service. Used CommonsShare web portal to interface with this,
using BDBags.

## Team Nitrogen -- FAIR data analysis and metadata descriptions
/Avi Ma'ayan/

Lots of different demos and prototypes: FAIRshake (scoring), BioJupies
(analysis), Datasets2Tools (metadata repository), FAIRsharing (repo for
standards), FAIRmission (entering medata). FAIRshake evaluation of data,
interoperable with different environments including general backend. Also
developed a game to learn the structure of the data commons.

## Team Oxygen -- Metadata mapping between models and search
/Anu Gururaj/

DATS4FAIR, a semantic metadta model supporting FAIRness. Other metadata model is
BioLinks: trying to map between them since they represent different parts of
datasets. Also have a pipeline for metadata ingestion and indexing: DATS2.2
format. Working on improving data search engine performance though metadata,
both phenotypic and genotyping data at the patient level. Also survey of
blockchain platforms and an explorer + immutable ledger.

## Team Phosphorus -- indexing and search of metadata
/Jonathan Crabtree/
Working on index and search of data, harmonizing metadata. Also helping with
organization and communication: slack and discourse for a helpdesk. Define
metadata definition and instantation: DATS. Nice overview of metadata model
instances and how they work together.

## Team Sodium -- organize data identifiers and retrieval
/Sarala Wimalaratne/
Working on Data Commons GUID requirements: identify fair digital objects.
Building GUID services: Identifiers.org DataCite, N2T: resolver, namespace and
object registration services. Full Stack teams use these services to get to the
data. GUID Core Metadata report: 3 types of identifiers Compate, ARKs and DOIs;
each one has a different set of core metadata. Namespace Service: register
namespaces for Compact Identifiers. Object Registration GUID services: generate
identifiers for individual datasets. Resolver Services across the multiple
existing services.

## Team Xenon -- FAIR2CURES full stack implementation
/Alison Leaf/

FAIR4CURES full stack platform built on SevenBridges, GUID Broker API. Focusing
on indexing and searching: level 1 -- metadata based, level 2 semantic
integration, level 3 -- aggregate data. For Level II, provides SevenBridges Data
Browser. Currently available on the Cancer Genomics Cloud, provides a nice way
to search, now moved over to FAIR2CURES platform.

## Team Hydrogen -- NIH data organization and management
/Chip Schwartz/

NIH, structuring how data looks and organized. Idea to get approaches into all of
the different stacks. Processes and tools for making this happen from a project
management suggestion.

## Data Stewards -- TopMed
/Josh Bis/
TopMed program identifying genetic variant associated with diseases, lots of
whole genome sequencing available on dbGAP, extensive phenotype and exposure
data. TopMed is a diverse set of samples: lots of phenotypes and different
ancestries. Build user narratives to center what the commons is. Provide tools
to find data and searching across phenotypic datasets, prepare genotypic
datasets including access, and tools for doing analysis.

## Data Stewards -- GTEx update
/Kristin Ardie/
Current public release available on GCP and AWS: V7 available now and V8 coming
out in summer of 2018. 17k RNA-seq and 838 WGS. V7 uses build 37, V8 moving to
build 38, uses a standard set of RNA-seq tools.

## Data Stewards -- GO and AGR
/Laurent-Phillipe Albou/
GO is a structured vocabulary for describing function of genes. AGR, facilitates
use of diverse model organisms. 20th anniversary of project. Working on new
pipeline to make query and usage easier. Nice use cases of GO enrichment APIs.

# Future talks

The second morning, groups presented their ideas for future work.

## Harmonizing the Data Commoners
/Rayna Harris/

A few of favorite things: linking disparate things, making analogies, and calls
to actions. Starts with discussion of talk from Avi at U Texas about Jupyter
notebooks and Harmonizome app. Her friend has a band and song/album called
Harmonizm so we get started with music thinking about Harmonization. Everyone
who is not here is missing an awesome photoshop of the album cover with Commons
people.

## Project Management Update
/Sasha Wait Zaranek/

Deliverables in April: crosscut data model proposal, demo on PIC-SURE API
and Gen3, harmonized APIs. Cross-consortium transparency on progress, clean
process for updating milestones, integration in deliverables and overall goals,
PMs have streamlined milestones for each group, demos and products evaluation
and sharing established (https://pilot.nihdatacommons.us/organize). Working out
improved lines of communication.

Really useful discussion around publicly sharing NIH Data Commons activities.
Generally want to do, but need to have a nice public facing message in addition
to the details and ongoing working steps.

## TOPMed
/Ben Heavner/

DCC -- Data Coordinating Center at TOPMed, coordinating informaticians and
scientists using the data. Reiterates 6 steps: developing analysis plan,
constructing phenotypic data set, preparing genomic data set, prepare variants
of interest, perform genotype-phenotype analysis, evaluate and refine results. A
lot of tricky decisions at each point which result in inconsistent
representation across studies. TOPMed whitepaper describing the work.

## GTex

Enabling search across Firecloud based on sample and donor annotations. Idea to
avoid transferring data; they have a lot. Need methods to expose manifests from
search so external users can make use of; want to use BDBag for this. Nice
discussion around how this will actually happen.

## Alliance of Genome Resources (AGR)
/Chris Tabone/
Model organism database focus: FlyBase, Mouse, Rat, Yeast, Worm, Zebrafish all
working together: http://aliancegenome.org. Data from all mods available in S3
in JSON format using a defined schema. Uses multiple ontologies (GO, SO...).
Stored in Neo4j, indexed in ElasticSearch, Federated data showed on AGR website.
Data Commons has ontologized, normalized, json. Schemas include Gene, Disease,
Orthology, Allele, Phenotype, Expression, Genotype. Want to figure out which
data does the commons want, should it be more cleaned up? Moving to BDBags for
data submission.

## TOPMed workflows registered on public site
/Cricket Sloan/
Provide sharable workflows for the consortium: aligner and variant caller from
UMich, hosted in GitHub. https://github.com/statgen/docker-alignment and
https://github.com/statgen/topmed_freeze3_calling Wrapped
into CWL and WDL workflows. Using a concordance tool in CWL to compare between
inputs. Dockstore is the public side for hosting them. Want to further develop
the checker tool, bigger test BED and exchange workflows.

## GUIDs and metadata

Next steps: core metadata for GUIDs coming in an RFC for comments. Represented
in schema.org but need to overlap with DATS model. Need to have a list of GUIDs
across all the different stacks and how to harmonize them.

## Working with harmonized and non-harmonized metadata
/Jonathan Crabtree/

Data Cube - 3 axes: Data Providers, Data Types, metdata object types. Have data
that is open/restricted plus harmonized/non-harmonized. Harmonization both
inside and outside of DATS model. Also working on replacing Data Portal using
the Crosscut Metadata model.

## FAIRshake redesign
/Daniel Clarke/

FAIRshake evaluates digital objects for FAIRness. Current code based uses
swagger to document and allow other tools to use API. Trying a different
approach going forward because changing the API is very difficult. New idea is
to put the API first, other tools can extend the API itself. The tool is then
modular and interoperable so can do different rubrics for fairness and ways to
create scores.

## Scientific Use Cases
/Jeremy Yang/

Data use case: using sex as a biological variable. Based around gene expression
data in GTEx and tissue specific expression. Demos allow exploration of
expression profiles. Second example using mouse data with alternative splicing.

## Data use cases with TOPmed WGS data
/Greg Versmee/

Working on TOPmed studies, 18 studies with 140,000 clinical phenotype variables
and 18,000 whole genomes. Plans: use of BDBags for integration, work on FAIR
principles. Collaborate through opening up the i2b2/tranSMART interface and
PIC-SURE API. Also have R tools for dealing with dbGAP. Next month will be
showing how to use the API in real life data.

## Proposal of first step interoperability demo with BDBags
/Carl Kesselman/

BDBags simple and doing something useful, very popular here at the Workshop.
Idea for another useful demonstration for interoperability: create BDBag with
simple data, ingest into full stack, compute on it, export BDBag and import into
a second stack. Demonstrate we can cycle data around between stacks and can talk
with everyone else. To get fancier, do an authenticated version.

## Connecting use cases to community and tool development
/Juergen Klenk/

Idea to connect design and development: collect use cases, derive personas,
design principles, develop APIs and tools. Really nice human centric approach to
turning use cases into code that actually solves the problem.

