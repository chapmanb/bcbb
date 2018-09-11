Collection of useful code related to biological analysis. Much of this is 
discussed with examples at [Blue collar bioinformatics][1].

All code, images and documents in this repository are freely available for all
uses. Code is available under the [MIT license](https://opensource.org/licenses/MIT)
and images, documentations and talks under the [Creative Commons No Rights
Reserved (CC0) license](https://creativecommons.org/share-your-work/public-domain/cc0/).

Some projects which may be especially interesting:

* CloudBioLinux -- An automated environment to install useful biological software and
  libraries. This is used to bootstrap blank machines, such as those you'd
  find on Cloud providers like Amazon, to ready to go analysis workstations.
  See the [CloudBioLinux][2] effort for more details. This project
  moved to its own repository at https://github.com/chapmanb/cloudbiolinux.
* gff -- A GFF parsing library in Python, aimed for inclusion into Biopython.
* nextgen -- A python toolkit providing best-practice pipelines for fully 
  automated high throughput sequencing analysis.  This project has 
  moved into its own repository: https://github.com/chapmanb/bcbio-nextgen
* distblast -- A distributed BLAST analysis running for identifying best hits in
  a wide variety of organisms for downstream phylogenetic analyses. The code
  is generalized to run on local multi-processor and distributed Hadoop
  clusters.

[1]: http://bcbio.wordpress.com
[2]: http://cloudbiolinux.org/
