Gad Getz -- Challenges in Cancer Genomics
Broad Institute medical and population genetics talk 26 Feb 2015

Gadi starts off with a description of what cancer is: initiating driver event,
passenger events, then more driver events that increase tumor
fitness. Differentiates clonal and sub-clonal events in the cancer population.
2 level approach to find cancer genes and pathways: characterization of
individuals with all changes, and population level interpretation. Power in
individuals is from depth of sequencing, in population is from number of
individuals. Existing cancer genome projects: TCGA/ICGC -- idea is to generate a
comprehensive catalog of cancer genes and pathways.

Characterization challenges to get complete base-level characterization of a tumor and
its evolutionary history:

- Flood of data. Need pipelines that scale and ways to store them. Firehose.
  Tools: muTect, Indelocator, SegSeq (CNVs) , dRanger and BreakPointer, PathSeq
  (pathogens). Most tools not publicly available:
  https://www.broadinstitute.org/cancer/cga/

- Mutation allelic fraction depends on purity, local copy number, multiplicity
  and cancer cell fraction. Need systematic benchmarking: target is false
  positives < 0.1/Mb. Two types of false positives: noise and germline events.
  Measure sensitivity and specifity using virtual tumors created from NA12878
  and NA12891 (parent). Okay, three types of false positives: also have
  cross-patient contamination: use ContEst to estimate.

- Different mutational patterns in different cancers. Example: AA > AC in
  esopageal cancer. Use NMF to segment samples into processes of
  mutations. Tumor types section by processes.

- Estimation of purity of cancer samples. Uses ABSOLUTE:
  https://www.broadinstitute.org/cancer/cga/absolute.

- Cancer cell fractions have large uncertainties. Cluster these to assign
  subclonal sections. Useful to look at this over multiple time points. At least
  two helps you assign evolution.

Second challenge is interpretation of populations. GISTIC -- scores regions
according to frequency and amplitude of copy-number events. MutSig -- score
genes according to number and types of mutations. Models background mutations
with a known distribution by chance. This assigns p-values to events per gene,
and can identify after correcting for multiple tests. Most genes mutated in a
small fraction of patients. Problem is that as sample size/mutation rates
increase, you get a larger number of genes that don't look right. Problem is
that we don't actually have an average background rate: different background
rate per gene. So heterogeneity at patient and base level need to take this into
account. There is also less mutation in highly expressed genes, and mutation
rate depends on gene synthesis time. MutSigCV takes these into account.

Uses 3 different signals to find cancer genes: number of mutations, clustering
of mutation, and in conserved sites. Did this analysis on both type specific and
overall across all cancers. Tumor Portal to explore these:
http://tumorportal.org Really nice.

Next question is if we've completed the catalog of tumor genes. No, still more
complex biology to detect. We're doing good on ones that are common (in 20% or
more patients) but genes in <2% of patients we're still at the
beginning. Estimation of sample power needed: 100,000 tumors (50 tumors x 2000
samples/tumor type). Currently have 15,000 cases.

Dark matter in cancer genomics. 48% lack a known driver mutation looking at
exome sequencing and CNV analysis, 38% lack a focal amplification in a driver
gene. Non-coding and epigenetic changes. Candidate cancer genes are still not
biologically verified with functional classification.
