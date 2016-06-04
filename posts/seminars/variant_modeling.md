## Variant modeling collaboration

I'm at a collaborative discussion at Broad about improving
[representation of variants](https://docs.google.com/document/d/1LX_1H1eHbkmisKps5WlFmE5FBJrhkVVEmxzAU0R2ySo/edit#heading=h.i017roj2lhe6).
The goals are to try and define how we manage and report variants to clinicians.
There are representatives from the following organizations:

- [FHIR Genomics](http://projects.iq.harvard.edu/fhirgenomics)
  _David Kreda, Gaston Fiore, Gil Alterovitz_
  FIHR provides clinical point of care data models
  through [SMART genomics](http://smarthealthit.org/). Helps with communication
  of clincal genomics data between EHR systems. FHIR is about data, and SMART is
  about creation of clinical genomics apps:
  [App Gallery](https://gallery.smarthealthit.org/). This helps manage the
  complication of [EL7 EHRs](http://www.hl7.org/EHR/). Clinical report ->
  clinical information -> sequences and variants. Sequencing resources then
  point out to full sequence sources like GA4GH.

- [ClinGen](http://www.clinicalgenome.org)
  _Bob Freimuth, Larry Babb_
  ClinGen provides a central authoritative resource for variants -- knowledge
  base. Genes, causitive variants, actionability, building genomic knowledge
  bases. Key goals around representing alleles and variant assessment. Supports
  ClinVar and many other sources. Practical output for this group:
  [ClinGen data model](http://datamodel.clinicalgenome.org/development/). Uses
  FHIR model for resources: Allele + Provenance. Focus areas: consistent variant
  representation, allele state, link knowledge evidence and interpretation to
  reference variant. Work to do: working definition of variaition, allele and friends.

- [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/)
  _Jennifer Lee_
  Variation resources at NCBI: ClinVar, dbSNP, dbVar, dbGap. NCBI 2.0 moving to
  AWS. Variation services: [Beacon](https://www.beacon-network.org/#/), Allele
  Registry and internal variation services. ClinVar data comes in HGVS or
  genomic coordinates: need an unambiguous allele registry to avoid confusion.
  Goal of variation service is to solve 95% of cases: machine readable
  nucleotide variations -- not SVs or unplaced/uncertain. Represent variations
  with simple method called FRV:
  Accession.Version:Position:Deleted-Sequences:Inserted-Sequence
  Computational definition -- everything removed and inserted. Workflow is HGVS
  to unambiguous FRV, then normalize as a canonical FRV. In the Allele registry
  have two ID spaces: contextual -- placement on transcript, canonical placement
  on reference genome.

- [GA4GH](http://ga4gh.org/)
  _Reece Hart_
  Reece talks about the goals of the GA4GH -- cooperating, global,
  decentralized. Framework is around federation. Some demonstration projects:
  Beacon, [Matchmaker Exchange](http://www.matchmakerexchange.org/). Reece has
  great [hgvs and uta](https://bitbucket.org/biocommons/hgvs) packages. He talks
  about the disconnect between analysis and reporting -- we don't push evidence
  from BAMs and VCFs downstream to clinical reporting. Horrible examples of
  things we fail on: different transcript representations of exons between UCSC
  and NCBI. Another tough example between genome and transcript -- repeats at
  splice sites so can represent variant in intron or exon, so will be in
  transcript or not depending on how you represent it.

The overall goals of the meeting are:

- Scope statement -- what variants do we cover? What are the useful outcomes?
- Good use cases -- submission, variant comparison
- Improving terminology -- ensure we're talking about the same things: alleles,
  variants. Synonyms, references, examples.
- Variant model: how do we represent variants? The ultimate goal -- start with
  representations from ClinVar, ClinGen and compare.
