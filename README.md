Collection of useful code related to biological analysis. Much of this is 
discussed with examples at [Blue collar bioinformatics][1].

Some projects which may be especially interesting:

* CloudBioLinux -- An automated environment to install useful biological software and
  libraries. This is used to bootstrap blank machines, such as those you'd
  find on Cloud providers like Amazon, to ready to go analysis workstations.
  See the [CloudBioLinux][2] effort for more details. This project
  moved to it's own repository at https://github.com/chapmanb/cloudbiolinux.
* gff -- A GFF parsing library in Python, aimed for inclusion into Biopython.
* nextgen -- Automated analysis pipeline for processing next generation
  sequencing data. This is tightly integrated with the Galaxy web framework.
* distblast -- A distributed BLAST analysis running for identifying best hits in
  a wide variety of organisms for downstream phylogenetic analyses. The code
  is generalized to run on local multi-processor and distributed Hadoop
  clusters.

[1]: http://bcbio.wordpress.com
[2]: http://cloudbiolinux.org/
