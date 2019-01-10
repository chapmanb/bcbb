## Sequencers

- Illumina NovaSeq, widely used

- PacBio Sequel, longer read lengths, 2-5ug of input required

- Oxford Nanopore: Promethion 70Gb/run (48 flowcells), Minion 15Gb/run max, in
  practice more like 1Gb unless you're really good (1 flowcell). Super long
  reads (>100kb wih record at 1Mb). Cool clever work: designed CRIPR every 400kb
  to cut genome. Oxford improved a lot in the last year. Some cool developments:
  direct sequencing of RNAs, generate cDNA for stability, 87% read identity,
  full transcripts, challenges with polyA tails but can estimate size in the
  future. Looked at [Lexogen SIRV spike in controls](https://www.lexogen.com/sirvs/)
  and good observed versus experimental.

## DeNovo assembly

Combo of two methods: 1. Long read input to scaffold on. 2. Use Hi-C like method
using Illumina; Dovetail genomics, [Phase Genomics](https://phasegenomics.com/hic-technology/)

## Long reads

Metagenomics -- used for assignment of 16S by linking to surrounding space.
Helps with assignment and function.

## Variant calling

Steve Lincoln at Invitae -- made synthetic variants off challenging human
variants to test protocols, then sent to 10 labs and highlighting variability in
reporting. Should be available now.

## Single Cell

sciSEQ from Schendure lab: does ATACseq, Hi-C, DNA, Methyl, RNA all as single
cell. Idea is to connect single cell chromatin architecture with disease
associated variations. If you have an association of variants with open
chromatin in a cell type then potentially functional in that cell type.

Rhapsody from BD, 19k cells with 3' DGE

Cell hashing/CITE-SEQ. Idea, add antibodies with tags onto cells and then
sequence through and get the tag. Can multiplex cells on the inputs; can also
detect duplicates with this and through them away. BD sells these antibodies and
platform agnostic. De-multiplexes multiple samples prepped together to reduce
prep costs.

10x two bead single cell method. Fuses two droplets (bead + bead and oligo).
Advantage is that you now that UIDs throughout the transcript instead of only at
3' end so can get at full transcripts. Improves access and costs for single
cell.

## Spatial genomics

Single cell loses the cellular positions, do not know the architecture from
where it came. Nanostring: slide decorated with oligo coupled to antibodies,
linked by light cleavable oligo. Slide stained by cell type and masted, then
release with light and examine. Vickovic at Broad uses microarray slide with
addressed oligos, 33x35 coordinates, did 1800 sections over the 33x35 to rebuild
the tissue.
