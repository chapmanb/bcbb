#!/usr/bin/env python
"""Install distblast software on every node of a Hadoop cluster.

This is an example of how to remotely add non-AMI data or software
to a Hadoop cluster kicked off with whirr.

Usage:
    cluster_install_distblast.py <cluster config file>
"""
import os
import sys
import subprocess

import fabric.api as fabric
import fabric.contrib.files as fabric_files

def main(cluster_config):
    user = "hadoop"
    addresses = _get_whirr_addresses(cluster_config)
    # software for all nodes on the cluster
    for addr in addresses:
        install_distblast(addr, user)
    # data on the head node
    dl_distblast_data(addresses[0], user)

def dl_distblast_data(addr, user):
    """Download distblast data from S3 bucket for analysis.
    """
    data_url = "http://chapmanb.s3.amazonaws.com/distblast.tar.gz"
    with fabric.settings(host_string="%s@%s" % (user, addr)):
        if not fabric_files.exists("distblast"):
            fabric.run("wget %s" % data_url)
            fabric.run("tar -xzvpf %s" % os.path.basename(data_url))

def install_distblast(addr, user):
    print "Installing on", addr
    with fabric.settings(host_string="%s@%s" % (user, addr)):
        work_dir = "install"
        if not fabric_files.exists(work_dir):
            fabric.run("mkdir %s" % work_dir)
        with fabric.cd(work_dir):
            distblast_dir = "bcbb/distblast"
            if not fabric_files.exists(distblast_dir):
                fabric.run("git clone git://github.com/chapmanb/bcbb.git")
                with fabric.cd(distblast_dir):
                    fabric.run("python2.6 setup.py build")
                    fabric.sudo("python2.6 setup.py install")

def _get_whirr_addresses(whirr_config):
    """Retrieve IP addresses of cluster machines from Whirr.
    """
    cl = ["/home/bchapman/install/java/whirr-trunk/bin/whirr", "list-cluster", "--config", whirr_config]
    return _addresses_from_cl(cl)

def _addresses_from_cl(cl):
    proc = subprocess.Popen(cl, stdout=subprocess.PIPE)
    stdout = proc.communicate()[0]
    addresses = []
    for line in stdout.split("\n"):
        parts = line.split()
        if len(parts) > 0:
            addresses.append(parts[2])
    return addresses

if __name__ == "__main__":
    main(*sys.argv[1:])
