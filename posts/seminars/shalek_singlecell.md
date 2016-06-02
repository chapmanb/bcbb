What can single cell genomics do for the clinic?
/Alex Shalek/

Motivation: genomes provide an average over a population of cells. Cellular
behavior is heterogeneous so this makes modeling over averages a problem.
Started working with single cell on TLR signaling with dentritic cells --
expected a single response but varies a lot between identical cells. Used
SMART-seq. Most variable genes are those involved in immune response,
constitutive stuff is the same. The useful thing is that genes co-vary across
single cells. In dentritic cells this corresponds to different maturity states.
Lots of information in this co-variation. Variability in sampling is small
compared to variability between single cells.

Glioblastoma -- 500 cells across 5 different tumor samples. First work was
trying to understand tumor heterogeneity. More heterogeneous tumors == poor
survival. See a large number of differences even in single gene like EGFR.

Tried to improve the methods by focusing on simpler systems, back to dentritic
cells. Want to dissect population-level responses -- there are some cells that
respond early (1 hour, instead of 4 hours). Early responders: precocious or
special? Special: you need communication to drive this response so the cells are
driving the expression of the rest of the cells. Emergent property from
interaction between cells.

Back to tumors -- look at entire system in a more systematic fashion. 19
patients, 5000 cells. Malignent cells from patients are very different, but
immune responses are similar. Use CNV calling from DNA, and then try to use RNA
expression to infer copy number: good correlation for large CNVs. Inferred RNA
results can distinguish malignant from non-malignant. Individual cells are
variable between high AXL and MITF, even though at the population level look
like a single responder. This explains clinical observations of resistance when
treated with single.

Similarity of immune system despite heterogeneity in tumors. CD2/3, CD4, CD8 ,
FOXP3, CD25 help differentiate different T cells. There is heterogeneity in
terms of which T cell exhaustion method is used. Environment around each tumor
in each patient is different -- no single approach despite apparent similarity.
With TCGA data, can look at cell-type specific genes to separate into groups:
T-cells versus CAFs. If you subset by tumor region, get different expression
across the space.

New approaches: cellular composition, distinguishing factors, microenvironment
and direct cell contacts. Cool idea of having a global picture of tissue
profiling. Need to define the cellular composition of complex tissue -- scale
single cell methods. Using droplets to increase throughputs. Drop-Seq: droplets
with barcoded beads. Every bead has PCR tag, cellular barcode and UMI tag. Every
oligo has a different UMI. To bring cells together device to combine beads and
barcodes and cells together. 95% specificity, 7000 genes, 12%. Evaluated on
retina cells to identify cell types: find the cell types we expect at resonable
ratios. Provides 39 clusters and new markers that you can validate.

DropSeq doesn't work for tumors -- not portable and complex for moving to
clinic. New approach is Seq-Well: stick things into a ton of tiny microwells.
Tests -- works very well and portable. Wells have positional information. Same
ratios from imaging and sequencing.
