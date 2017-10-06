I'm at the [DNAnexus connect](http://www.dnanexusconnect.com/) user meeting in
Cambridge. These are my notes from Oct 5th, the first day of the conference.
Richard Daly starts the day off with a high level overview of where DNAnexus is
planning to go over the next year. He emphasizes the importance of workflow
flexibility, data organization, sharing and management.

## Pharma 2020 -- a look ahead in a changing world
/Steve Arlington, Pistoia Alliance/

Goal is to provide a vision of the future and Steve describes the difficulty in
predicting ahead based on the past. Predicting that the world medicine will
double over time, but not clear where and how that will happen. A challenge is
that existing medicines are good and need to provide better, cheaper and more
effective treatments. Drivers of additional needs: aging population, obesity,
cardiovascular disease. More people live longer with chronic diseases.`HONDAs
(Hypertensive, Obese, Non-compliant, Diabetic Asthmatics) generate 70% of
healthcare costs. Challenge is sustaining costs; pressure from regulators to
have ways to do this. Need to prove that you're adding value. Currently patients
don't have recourse if approach/drug doesn't work and we should be able to do
better. Research and development productivity has declined over the past 15
years and costs per approved molecule are too high. Healthcare continues to
increase in cost due to both this and larger number of people to treat.
Approach: provide a new lever by focusing directly on patient outcomes.
Understanding and incorporating health economics and phenotypic data.
Integrate the full pipeline you'll need for approval at an early stage. Need to
incorporate technical advances: lots of unstructured data that needs
integration. Example: biosensor that records everything you do and calorie
consumption. Could do now, but when and where are we going allow it? How can we
integrate it and improve behavior? Shift is moving from doing everything
yourself to integrated collaborative companies.

## Moving from a Reactive to a Proactive Approach to Cancer Care
/Bill Dalton, M2Gen/

[ORIEN (Oncology Research Information Exchange Network)]
(http://www.oriencancer.org/): goal is to bring together multiple cancer centers
to share and collaborate. Principles: inclusiveness, accessibility,
public-private partnerships and collaborative. There is a lot of overlap with
work in the open source community on
[FAIR data principles](https://www.force11.org/group/fairgroup/fairprinciples).
Total Cancer Care approach: discover through access to molecular data and
collaborative data; development by matching patients with trials, and delivery
by connecting pharmas and patients. The overall goal is to improve patient
orientation. Consenting approach moved from face-to-face to a 6 minute video.
Approved by many IRBs in multiple states. In terms of usage, many are happy to donate
clinical data for research usage. A technical challenge is storing and managing
the data. M2Gen does the work of making data available via warehousing. Approach
is to collaborate and publish on outcomes to make it more available. Trying to
use data to anticipate needs up front and know options up front. Working on a
single IRB and trial agreement to make it easy to act on trial options: want to
have both best options and ability to execute on it. Platform is hosted on AWS
with molecular hosting via DNAnexus -- then exposed as a data warehouse.
Focusing on immunooncology using signatures based on things like mutation load
to assign patients. Summary: move from reactive care to proactive care. ORIEN
has a patient advocacy group driving goals. using signatures based on things
like mutation load to assign patients. Summary: move from reactive care to
proactive care. ORIEN has a patient advocacy group driving goals.
[Patient facing website for Total Cancer Care](http://orientcc.org/). Biggest
current challenge: getting and harmonizing clinical data.


## Birth of the Cool: Miles Past 200K Exomes
/Jeff Reid, Regeneron/

[Regeneron Genetics Center](https://www.regeneron.com/genetics-center) focus is
to create an environment where studies interoperate: families, big population
studies, phenotype specific cohorts. Collaboration with UK Biobank -- 500k
participants and brilliant open resources moving from genotyping to sequencing.
Idea it to apply automation to the entire process of sequencing and analysis:
100% cloud on top on DNAnexus. Highlights recent work
[sequencing 100k samples in Geisinger DisocvEHR](https://www.biorxiv.org/content/early/2017/10/03/197889);
unique finding is that there is a higher degree of relatedness than expected.
People live together so makes sense but not taken advantage of or considered in
large scale studies. Can use relatedness to phase genomes and understand
compound heterozygous knockouts. The theme here is fully exploiting it, you now
can increase a larger percentage of knockouts due to the phased hets. In 200k
cohort, 23% of genes have a homozgous knockout and will get higher with compound
het knockouts. Do all-by-all association results for screening; pre-computed so
you can lookup any association in a big database.

Highlights some cool work that can do with this data. Scan for detectable
chromosomal abnormalities based on allele balance anomalies. Example, see
alleles at 0%, 25%, 50% indicating there is something non-diploid. Nice parallel
to depth based analyses using AF as well. Can see runs of homozygosity, sub
chromosomal level deletions. There is a big diversity of chromosome X
anomalies and see a lot of 1:1000+ events that patients don't know about.

Find clinically actionable findings in 1:29 participants. Working on ways to
feed this back to patients. Example of Barbara Barnes -- found BRCA1 mutation
that completely changed treatment. One challenge is figuring out IRBs to allow
reporting back for CNV and other issues.

## Integrative Informatics: How to Bring Semantic and Tabular Data Together to Drive Insight
/Mathew Woodwark, MedImmune/

Matthew talks about data integration work ongoing at MedImmune and AstraZeneca.
Overall goal is to provide targeted therapies. Integrating across a whole ton of
different input data sources. Focus is on answering scientific questions. Avoid
enterprise data warehouses. Goal is to better integrate genetics throughout the
pipelines; want the underlying genetic reason for the change to help stratify
patients. Emphasis on paying for success so understanding mechanism is key.
Needs: better ways to represent population level data -- compute to data and
integrates results. Integrative informatics: it can be integrated but only bring
it together when you want to answer the question. Change in goals and score:
answer questions, not integrate everything. Changes the focus and what you want
to do. Underlying tech is a semantic graph, allowing deconstruction of
biological questions to graph. Advantages of linked data: integrating with other
data providers that are already doing this like EBI.

## Enabling Standard Bioinformatics Processing at UPMC’s Genome Cente
/Annerose Berndt, University of Pittsburgh Medical/

[University of Pittsburgh Medical (UPMC)](http://www.upmc.com/Pages/default.aspx)
huge integrated health care provider with fully integrated system in multiple
countries. The [UPMC Genome Center](http://path.upmc.edu/genome/Index.htm) is a
year old attempt to combine the disparate genomic work happening across the
institution. Common theme:
coordinated consent across biobanks and multiple participants. Other challenges:
multiple EMR systems that feed into LIMS
They wanted to give a central repository for genomic data but avoid
having an internal footprint, hence use of DNAnexus on AWS. Stream data from
sequencers directly to AWS, so have direct connect. Why did they decide on
DNAnexus? Security and compliance, analysis automation using Sentieon and Edico.
Timing: 12 hours bcl to fastq, 2 hours fastq to bcl. Two focuses now: cancer and
pediatrics.
