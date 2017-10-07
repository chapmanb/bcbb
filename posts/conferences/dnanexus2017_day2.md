I'm at the [DNAnexus connect](http://www.dnanexusconnect.com/) user meeting in
Cambridge. These are my notes from the second day, Friday Oct 6 2017 building
off the [notes from day 1](https://github.com/chapmanb/bcbb/blob/master/posts/conferences/dnanexus2017_day1.md).

## State of Precision Medicine: Update
_David Shaywitz, DNAnexus_

David sets the context of talks for the day with history of genome hype: knowing
the genome, NGS + cloud + multiple genomes. Skepticism needed because
progress is hard, and driven by champions that persist in moving forward on the
right problems. Core mission of DNAnexus is powering champions of genomic medicine.

## Next Generation Science: Developments and Roadmap for DNAnexus Science
_Andrew Carroll, DNAnexus_

Andrew leads the science team at DNAnexus talking about work from the past year
and plans for the upcoming year. Promises from last year: DNAnexus easier to
use, improve job monitoring, improve user management, expand and make tools more
widely available, collaborations with external scientists. Tools: add ~1 app per
week over the last year. Examples: Sentieon, Edico. Easier to use: added support
for [Common Workflow Language](https://github.com/dnanexus/dx-cwl) and
[Workflow Definition Language](https://github.com/dnanexus-rnd/dxWDL). Also
support containers with
[dx-docker](https://wiki.dnanexus.com/Developer-Tutorials/Using-Docker-Images).
For in depth debugging work can use
[workstations](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations)
to get and snapshot in progress work. New feature: `dx notebook` which launches
web notebooks from the commandline with snapshots. Great improvements to
visualization from [Maria Nattestad](https://github.com/MariaNattestad).
DNAnexus joint caller is [GLnexus](https://github.com/dnanexus-rnd/GLnexus) --
working in practice for real samples.

Plans for the next year: improve communication of new tools and features, focus
on the app development process to make it easier, expand scientific
collaborators with customers and the community.

## Microsoft’s Vision: Improving Global Health through Digital Transformation
_Clifford Goldsmith, Microsoft Health and Life Sciences_

Microsoft has an impressive global mission to help, focusing in this talk on the
US Health strategy. Better health, better care, productivity and lowering costs.
Practically improving the clinical experience. Focusing on genomics as part of a
medical patient journey. Useful way to frame work around patients directly. With
move from paper to digital can personalize the journey for patients:
crowdsourced version of implementation.

## Cancer Genomic Analysis Tools and Pediatric Genome Data in the Cloud
_Jinghui Zhang, St. Jude Children’s Research Hospital_

Pediatric cancer challenges: unique properties relative to adult cancer. Need to
make data available for building tools: 700 paired tumor/normal WGS, 1200
exomes + RNA-seq. Data available through the
[Pediatric Cancer Genome Project (PCGP)](https://www.stjude.org/research/pediatric-cancer-genome-project.html).
Key findings: only 45% of driver genes overlap with adult cancers. Structural
changes account for 62% of drivers -- importance of WGS. Great customized
visualization tools: [ProteinPaint](https://www.stjude.org/research/shared-resources/technology-licensing/technologies/proteinpaint-web-application-for-visualizing-genomic-data-sj-15-0021.html).
Integrate mutation, gene fusion and expression. Challenge for sharing data with
the global community: raw data in EGA, code/pipeline sharing: GitHub, cloud
distribution by machine learning. Emphasizes the support burden. Goal with
moving to the DNAnexus wsa to provide a platform for running tools and making
data available. Challenge in transitioning from research to clinical genomics:
specific FFPE pipelines. For detection of fusions: CICERO, coupled with manual
inspection using visualization tools. Amazing examples of complex visualizations
that are difficult to detect without de-novo assembly. St Jude Cloud is an
ecosystem of cloud-based offerings, making use of existing resources on
DNAnexus like the [PeCan data portal](https://pecan.stjude.org/). Contains data
from a worldwide set of collaborators.

## Introducing DNAnexus Clinical Trial Analytics and Data Management PaaS
_Omar Serang_

Current state of clinical development: more data but more to look through,
trials are incredibly expensive and a fragmented set of solutions. Use of NGS in
clinical trials is on the rise. Increased success rate for oncology drugs when
using pre-existing biomarkers to stratify. DNAnexus providing clinical trial
data services. The amount of regulatory requirements for clinical trials are
incredible, nice to have someone taking care of this and ensuring the platform
is compliant. Stresses importance of validation highlighting precisionFDA.

## The Janssen Clinical NGS Platform
_Ling Yang Hao, Janssen Pharmaceutical_

Focus is on drugs for Inflamatory Bowel Disease (IBD): Ulcerative Colitis and
Crohn's disease. Due to chronic inflammation of the GI tract. Looked at genetic
component and different between different populations: European versus Asian.
Different risk alleles and complex genetic association -- not exactly sure how.
Promising ideas: microbiome based approach to improve collection of good
bacteria relative to bad. VE202 is the drug product to approach this problem but
not clear how to analyze microbiome data as a drug. Implemented analysis
pipeline for stool samples on DNAnexus with configurable testing of different
parameters.
