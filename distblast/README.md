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

2. Install whirr (currently you need to [build whirr from svn][1a];
   version 0.3.0 is missing some patches to allow use of custom AMIs):

    sudo apt-get install whirr

3. Create ~/.ec2/distblast.whirr describing connection information.
   Use a Amazon machine image (AMI) with the necessary software to
   run your application; for instance CloudBioLinux.

    whirr.service-name=hadoop
    whirr.cluster-name=distblast
    whirr.instance-templates=1 hadoop-jobtracker+hadoop-namenode,1 hadoop-datanode+hadoop-tasktracker
    whirr.provider=aws-ec2
    whirr.identity=${env:AWS_ACCESS_KEY_ID}
    whirr.credential=${env:AWS_SECRET_ACCESS_KEY}
    whirr.private-key-file=${sys:user.home}/.ssh/id_rsa
    whirr.public-key-file=${sys:user.home}/.ssh/id_rsa.pub
    whirr.hardware-id=t1.micro
    whirr.location-id=us-east-1c
    whirr.hadoop-install-function=install_cdh_hadoop
    whirr.hadoop-configure-function=configure_cdh_hadoop
    whirr.login-user=ubuntu
    whirr.image-id=us-east-1/ami-4e57a227
    jclouds.ec2.ami-owners=678711657553

4. Start up the cluster and login:

    whirr launch-cluster --config ~/.ec2/distblast.whirr

5. Install distblast on each node and organism data on the head node:

    python2.6 bcbb/distblast/ec2/cluster_install_distblast.py ~/.ec2/distblast.whirr

6. Run a distributed BLAST from your local machine

 a. Start the proxy in a separate terminal. Allows you to connect to
 the cluster and also view the Hadoop web interface in
 http://cluster_machine_name:50030/. See the proxy information in the
 [quick start guide][4] for details on setting up your browser to use
 the proxy.

    % . ~/.whirr/distblast/hadoop-proxy.sh

 b. Setup enviornment to use local Hadoop pointed to remote cluster

    % export HADOOP_HOME=/usr/lib/hadoop
    % export HADOOP_CONF_DIR=~/.whirr/distblast
    % wget http://chapmanb.s3.amazonaws.com/distblast.tar.gz
    % tar -xzvpf distblast.tar.gz
    % cd distblast
    % python2.6 bcbb/distblast/hadoop/hadoop_run.py \
      bcbb/distblast/hadoop/distblast_streaming.py \
      org_configs/test.yaml base_config.yaml input output

7. Finished: logout, terminate the nodes and remove the cluster:

    whirr destroy-cluster --config ~/.hadoop-cloud/distblast.properties

[1]: https://wiki.cloudera.com/display/DOC/Whirr+Installation
[1a]: https://cwiki.apache.org/confluence/display/WHIRR/How+To+Contribute
[3]: https://wiki.cloudera.com/display/DOC/Hadoop+Installation+(CDH3)
[4]: http://incubator.apache.org/whirr/quick-start-guide.html

### Local hadoop cluster on ubuntu

Follow [installation documentation][2]:

1. Start hadoop

    sudo -u hadoop bash
    /usr/lib/hadoop/bin/start-all.sh

2. Run scripts

    % cd distblast_data
    % python2.6 bcbb/distblast/hadoop/hadoop_run.py \
      bcbb/distblast/hadoop/distblast_pipes.py \
      org_configs/test.yaml base_config.yaml input output

#### Debugging tips

Default timeouts are set to 10 minutes, which will result in Hadoop
sitting for long periods of time when something is wrong. While
debugging issues, you can reset this with `mapred.task.timeout` in
`mapred-site.xml` to get informative error messages quicker.

The MapReduce job tracker web interface at http://localhost:50030/ has
lots of useful error information you may not be seeing from the
commandline.

[2]: http://www.michael-noll.com/wiki/Running_Hadoop_On_Ubuntu_Linux_(Single-Node_Cluster)
