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

Follow [Cloudera script documentation][1]:

1. Install Cloudera's Hadoop distribution, following [Cloudera's instructions][3].

2. Install whirr:

	sudo apt-get install whirr

3. Create ~/.hadoop-cloud/distblast.properties describing connection information.
   Use a Amazon machine image (AMI) with the necessary software to
   run your application; for instance CloudBioLinux.

	whirr.service-name=hadoop
	whirr.cluster-name=distblast
	whirr.instance-templates=1 jt+nn,1 dn+tt
	whirr.provider=ec2
	whirr.identity=yyy
	whirr.credential=yyy
	whirr.private-key-file=/home/chapmanb/.ec2/id-sobchak.keypair
	whirr.hardware-id=t1.micro
	whirr.location-id=us-east-1c
	whirr.hadoop-install-runurl=cloudera/cdh/install
	whirr.image-id=us-east-1/ami-4e57a227
	jclouds.ec2.ami-owners=678711657553

4. Start up the cluster and login:

	whirr launch-cluster --config ~/.hadoop-cloud/distblast.properties

5. Install distblast on each node:

	python2.6 bcbb/distblast/ec2/cluster_install_distblast.py  ~/.hadoop-cloud/distblast.properties

6. In AWS console (https://console.aws.amazon.com/ec2/):

         - Create data volume from snapshot
         - Attach data volume to head node as /dev/sdf

7. Login to the cluster and mount the data volume containing the organism data:

	whirr list-cluster --config ~/.hadoop-cloud/distblast.properties
	ssh -i ~/.ec2/id-sobchak.keypair ubuntu@first-ip-address
        # sudo mkdir /mnt/phyloblast
        # sudo mount /dev/sdf /mnt/distblast

8. Run the cluster

	# cd /mnt/distblast/
        # python2.6 ~/install/bcbb/distblast/hadoop/hadoop_run.py ~/install/bcbb/distblast/hadoop/fasta_process.py \
           org_configs/test.yaml base_config.yaml input output

[1]: https://wiki.cloudera.com/display/DOC/Whirr+Installation
[3]: https://wiki.cloudera.com/display/DOC/Hadoop+Installation+(CDH3)

9. Finished: logout, terminate the nodes and remove the cluster:

	whirr destroy-cluster --config ~/.hadoop-cloud/distblast.properties

#### Old way -- using deprecated python script

1. Install from whirr-trunk/contrib/python

2. Create ~/.hadoop-cloud/clusters.cfg with configuration called distblast

3. Start up cluster:

	hadoop-ec2 launch-cluster --user-data-file bcbb/distblast/ec2/hadoop-ec2-init-remote.sh \
			          distblast 1 nn,snn,jt 1 dn,tt
	hadoop-ec2 login distblast
	ssh -i ~/.ec2/id-sobchak.keypair root@ec2-50-16-13-181.compute-1.amazonaws.com

4. Finished, remove cluster:

        hadoop-ec2 terminate-cluster distblast
        hadoop-ec2 delete-cluster distblast

### Local hadoop cluster on ubuntu

Follow [installation documentation][2]:

1. Start hadoop

         sudo -u hadoop bash
        /usr/lib/hadoop/bin/start-all.sh

2. Running pydoop

        export HADOOP_HOME=/usr/lib/hadoop

[2]: http://www.michael-noll.com/wiki/Running_Hadoop_On_Ubuntu_Linux_(Single-Node_Cluster)
