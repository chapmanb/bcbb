## Overview

Example distributed BLAST workflow that can run in parallel on a single
multi-processor machine or fully distributed on a Hadoop cluster. The
application is to compare a reference genome in FASTA format against many
different organism databases. Using BLAST, the best sequence similarity 
hit in each database is identified and summarized in a tab separated output
file.

The job can be parallelized by input organism, by each protein in the organism
genome, and by each BLAST comparison. This provides an interesting example
application for demonstrating Hadoop and multi-processor parallelization. It
also lends itself to cloud computation on Amazon EC2 resources, allowing 
distribution as a ready to replicate experiment with source controlled 
scripts and publicly available data volumes.

## Working notes

### Initializing a Hadoop cluster with EC2

Follow [Cloudera script documentation]:

1. Install script hadoop-ec2

   wget http://archive.cloudera.com/cdh/3/hadoop-0.20.2+320.tar.gz
   cd hadoop-0.20.2+320/src/contrib/cloud/src/py
   python setup.py install

2. Create ~/.hadoop-cloud/clusters.cfg describing connection information.
   Use a Amazon machine image (AMI) with the necessary software to
   run your application; for instance CloudBioLinux.

3. Start up the cluster and login:
   hadoop-ec2 launch-cluster small-cluster 1
   hadoop-ec2 login small-cluster

4. In AWS console (https://console.aws.amazon.com/ec2/):
     - Create data volume from snapshot
     - Attach data volume to head node as /dev/sdf
     In terminal on head node:
     - Get data
     # mkdir /mnt/phyloblast
     # mount /dev/sdf /mnt/phyloblast
     - Install required scripts
     # git clone git://github.com/chapmanb/bcbb.git
     - XXX need to do this on every node or put in image:
       # cd bcbb/phylogenetics/blast/
       # python2.6 setup.py install
     - Run:
     # cd /mnt/phyloblast/
     # python2.6 ~/bcbb/hadoop/hadoop_run.py ~/bcbb/hadoop/fasta_process.py org_configs/test.yaml base_config.yaml input output

[1]: https://wiki.cloudera.com/display/DOC/CDH+Cloud+Scripts

5. Finished, terminate the nodes and remove the cluster:
   hadoop-ec2 terminate-cluster small-cluster
   hadoop-ec2 delete-cluster small-cluster

### Local hadoop cluster on ubuntu

Follow [installation documentation][2]:

1. Start hadoop
   sudo -u hadoop bash
   /usr/lib/hadoop/bin/start-all.sh

2. Running pydoop:
   `export HADOOP_HOME=/usr/lib/hadoop`

[2]: http://www.michael-noll.com/wiki/Running_Hadoop_On_Ubuntu_Linux_(Single-Node_Cluster)
