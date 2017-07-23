I'm at the
[2017 Bioinformatics Open Source Conference (BOSC)](http://www.open-bio.org/wiki/BOSC_2017)
in Prague. These are my notes from the day 2 morning sessions, including the
epigenetic personalized medicine keynote, community building, citizen science
and late breaking talks.

BOSC is a two day annual conference organized by the
[Open Bioinformatics Foundation](https://www.open-bio.org/wiki/Main_Page). It's
an open community interested in building reusable and interoperable tools for
doing biological research.
All video from these presentations will be available from the
[schedule page](http://www.open-bio.org/wiki/BOSC_2017_Schedule)
and slides from the [f1000 channel](https://f1000research.com/collections/bosc).

# ISMB keynote

## Bioinformatics for Personalized Medicine: Looking Beyond the Genome
_Christoph Bock_

[Talk focuses on how the epigenome will impact personalized medicine](http://medical-epigenomics.org/meg/).
Most drugs only help a fraction of people, because each patient has an
individualized molecular background. In cancer, standard of care right now
includes gene panel or targeted sequencing to match patients to therapies. Many
more challenges ahead to integrate multiple data types, model heterogeneity and
quantify and propagate prediction confidence.

Research focus is epigenomic alterations: methylation, chromatin marks. Builds
off cancer genetics looking at driver and passenger mutations. At the next
level, inherited epigenomic reprogramming including complex interactions with
the genetics.

Focus is on Ewan sarcoma, an aggressive bone tumor in children. It is
genetically simple cancer, driven by a single EWS-ETS gene fusion event with few
additional alterations. Hypothesis is that epigenetic heterogeneity may explain
clinical response heterogeneity. Used bisulfite sequencing. As a first pass,
unsupervised methods found no useful clinical relationships. Nathan Sheffield
dug further into the data, [developing a method for finding regulatory footprints
from methylation](www.nature.com/nm/journal/v23/n3/full/nm.4273.html). Outcome:
has a epigenetic signature of initial cell type.

Second example, cancers of unknown primary site. Profile tumor samples and try
to use methylation marks to predict the initial site of the tumor. This has
validated clinical utility, taking 5 years to reproduce and translate to clinic.
Next step is a clinical trial, but promising step for better assigning
therapies. So, cells contain an epigenetic record of development history and can
use to improve treatments.

Next steps are to try and define an epigenetic landscape over time.
[International Human Epigenome Consortium (IHEC)](http://ihec-epigenomes.org/):
1000 reference genomes covering cells over entire human body. Can we use this to
infer differentiation hierarchies from epigenetic data. Can infer computational
lineage tree that matches FACS sorted manual work demonstrated using blood
differentiation. [LOLA: bioconductor package](https://github.com/nsheff/LOLA)
that provides GSEA like method for genomic region BED files, used to
differentiate cell types. Lots of great tool development from the lab.

Second example: using these methods to stratify leukemia patients -- chromatin
signatures associated with different disease subtypes. Chromatin data gives
principal components capturing the variability in the disease, placing patients
on a landscape. Idea of a progression zone, how they change over time.
Incredibly useful analogs to genetic heterogeneity analysis. Cancer, and biology
in general, is staggeringly complex.

Next idea, how do drugs modify the chromatin/epigenetic landscape. Try to model
using CRISPR for high throughput epigenomic editing. Problem: need to introduce
1000s of alterations into cancer cells. Approach: CRISPR single cell screening:
CRISPR library, infect in bulk, apply drug/virus stimulus, perform single cell
RNA-seq, link transcriptome to CRISPR alteration. Method is
[CROP-seq](https://github.com/epigen/crop-seq/). Similar to
[Perturb-Seq](https://en.wikipedia.org/wiki/Perturb-seq) which uses an expressed
barcode to tag the gRNA. CROP-seq meant to scale beyond thousands of samples.
Brilliantly open work with protocols available online.

Summary: epigenomic marks are generalization of cell type, capturing development
history and future potential of cells. Prediction methods to understand require
embracing and understading complexity. Practical clinical work: keeping patients
out of danger zone of disease progression. Need to engage the patients: this is
a cornerstone for personalized medicine. So nice to see discussion of citizen
science, open discussion and focus on education.

## Community Building and Citizen Science

## BeerDeCoded: exploring the beer metagenome
_Jonathan Sobel_

[BeedDeCoded](http://www.genome.beer/) is a personal project aimed at a general
audience. Project at the cool [Hackuarium](http://www.hackuarium.ch/en/)
hackerspace in Lausanne. Idea is to work with drinkers, homebrewers and
microbreweries. Trying to engage at festivals and events -- it's a way to engage
with the public and discuss issues around DNA sequencing. Crowdfunded on
Kickstarter to fund the project: nominate a beer, sample it, use
[BentoLab](https://www.bento.bio/bento-lab/) to prepare and sequence with
Nanopore MinION. Did beer extraction workshops. Practically, aligned against ITS
database. All
[analysis and code on GitHub](https://github.com/beerdecoded/Beer_ITS_analysis).
Results: found higher variety of specifies in complex fermented beers.
Visualized with cluster plot of beers and organisms. Nice example of making
science accessible to the public in an interesting way.

## Supporting curation communities & collecting technical dividends
_Monica Munoz-Torres_

[Apollo](http://apollo.berkeleybop.org/) is a great curation tool for improving
genome representations. It's a multi-user, graphical web application for edit
and improve annotations. It's also a social network for curators. On the
architecture side, uses jBrowse and then communicates with a backend Apollo
server. It now has a web services API that interacts with Galaxy. Idea is to
have more genomes with better annotations and have a great user base using
Apollo for this. Curation improves when there is dialog between and within
communities. Tons of great examples of different genome communities working with
Apollo, from Koalas to Heliconius butterfies to Arthropods. Communities are both
users and contributing back to Apollo code base.

## Journal of Open Source Software (JOSS)
_Pjotr Prins_

The [Journal of Open Source Software (JOSS)](http://joss.theoj.org/) initiated
based on observation that software developers need a track record of
publication. However most software developers don't focus on this since it ends
up being less of a priority next to coding, developing and documenting. Code is
often a better demonstration of a project than a paper. JOSS is all on GitHub
and papers are short descriptions of the software. JOSS review process is
lightweight, fast, free and public. It's fully indexed and a real journal.
Important point made during Q/A: DOIs don't mean anything with regards to being
peer reviewed.

## Building an open, collaborative, online infrastructure for bioinformatics training
_Bérénice Batut_

There is an increasing demand for learning bioinformatics. Currently happens
through informal training but surveys show most would prefer online or local
training options. [Galaxy](https://usegalaxy.org) is a nice scalable platform to
get folks to learn bioinformatics.
[Galaxy Training Network](https://galaxyproject.org/teach/gtn/) aggregates
training resources, more than a registered event per month over the past few
years so lots of ongoing training with this work worldwide. Difficulties: hard
to keep up with improvements and changes to Galaxy and tools. Established a
shared online set of tutorials. Very nice markdown to tutorial HTML process.
Also collecting easy to use small datasets for tutorials. Includes required
tools in Galaxy needed for each tutorial to make it easier to install and setup.
They've really thought through all of the shared pain of building up training
materials. Galaxy training materials all available from
[TESS](https://tess.elixir-europe.org/) and tagged with metadata. Incredible
resource to help with scaling teaching to meet the needs of researchers.

## Software and social strategies for community sourced biological networks and ontologies
_Dexter Pratt_

Aim is to get more and better exper content by engaging with the community.
[NDEx](http://www.ndexbio.org/) is infrastructure for finding ontologies and
other networks. Idea is to engage lab researchers to keep things up to date.
Not just curators to have user driven content with discussion. Idea is to create
topic focused communities using social strategies to motivate: making it easy to
publish by providing collections of networks. Provide full provenance and adding
automation to improve networks and ease burden on researchers.

## Distance-based, online bioinformatics training in Africa: the H3ABioNet experience
_Nicola Mulder_

[H3ABioNet](http://www.h3abionet.org/) builds infratructure for Africa, talking about fulfilling the need
for basic bioinformatics training for scientists. Provide a distributed
classrooms models, with pre-recorded lectures and practical assignments. All
material is creative commons and freely available. Used in a variety of
classrooms across Africa -- interact within classes and also between locations.
Provides virtual classrooms and live sessions. Amazing interactivity and
scalability.

# Late-Breaking Lightning Talks

## Recent object formation in the core of Galaxy
_Martin Cech_

An update on the work in the [Galaxy](https://usegalaxy.org) community in the
past year. Martin shows an amazing list of public Galaxy servers, focuses on
internationalization work. Emphasizes the reproducibility stack with
bioconda and biocontainers: deprecated the Galaxy package manager in favor of
bioconda. Jupyter based notebooks embedded inside Galaxy. Other interactive
tools: Phinch, RStudio, Neo4j. An incredible ecosystem of diverse users.

## Reproducibility of computational workflows is automated using continuous analysis
_Brett Beaulieu-Jones_

Nice followup on a preliminary talk at BOSC 2017, focusing on how to improve
reproducibility of workflows. Big issue, half of papers don't cite version of
software used.
[Continuous analysis](https://github.com/greenelab/continuous_analysis) provides
a continuous integration based approach to running analyses. Provides
reproducibility and can see changes of analysis over time. Can share audit logs
even if can't share original data. Great analysis with Geisinger health on 750k
patients how a code change influenced results.

## Full-stack genomics pipelining with GATK4 + WDL + Cromwell
_Kate Voss_

[GATK](https://github.com/broadinstitute/gatk) is a widely used variant calling
toolkit from Broad that is now open source as of GATK4. Yay.
[WDL](https://software.broadinstitute.org/wdl/) describes workflows using a DSL.
[Cromwell](https://github.com/broadinstitute/cromwell) is the workflow runner
that uses WDL. GATK4 is 5x faster than GATK3 and incorporates Spark. Now publish
best practice practices for somatic variations, and germline/somatic copy number
changes. Broad sequencings 24Tb of data a day. Cromwell runs on Google Compute
Engine and working on AWS + Alibaba.

## ToolDog - generating tool descriptors from the ELIXIR tool registry
_Kenzo-Hugo Hillion_

[bio.tools](https://bio.tools/) is a registry for finding tools which helps with
visibility and discoverability. [ToolDog](https://github.com/bio-tools/ToolDog)
uses a unique source of information in bio.tools to generate things like EDAM
and Galaxy datatypes. ToolDog inspects the source code using Docker generating a
template for a tool description. It can use existing descriptions and annotate
them. Great way to improve tools with descriptions.

## BioThings SDK: a toolkit for building high-performance data APIs in biology
_Chunlei Wu_

[BioThings SDK](http://docs.biothings.io/en/latest/) builds off BioThings talk
yesterday, describing the API used for building on top of it. Backend of
BioThings is a DataHub, API, then front end. Two ways to use BioThings SDK.
First, turn a data file into a BioThings API people can use: parse into JSON,
describe how to index, then serve out as JSON with Elasticsearch. Second,
aggregate multiple data sources and keep them up to date: parse multiple inputs,
use MongoDB for staging.

## Integrating cloud storage providers for genomic analyses
_Ted Liefeld_

[GenomeSpace](http://www.genomespace.org/) makes computational analyses
integrated by connecting tools and data. 20 different tools integrated and
multiple places to store files: S3, Dropbox, Google Drive, SwiftStack. Moves
data around between different locations. Provides a FileSystemInterface Java SDK
to access custom providers if want to use a non-provided storage solution.
Provide a REST-ful API for retrieval.

## Fighting Superbugs with Open Source Software
_Kai Blin_

Superbugs are drug resistant microbes, shows depressing and scary pictures of
rate of drug resistance across European countries. These occur because of
continued misuse of antibiotics. Use bacterial natural products which produce
useful antibiotics and are very hard to synthesize. Nice to see some chemical
structures. Idea: can we use genomes to identify these natural products. Open
source software to do this: [anitSMASH](https://bitbucket.org/antismash/antismash).
