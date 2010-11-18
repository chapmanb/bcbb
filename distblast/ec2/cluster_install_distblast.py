#!/usr/bin/env python
"""Install distblast software on every node of a Hadoop cluster.

This is an example of how to remotely add non-AMI data or software
to a Hadoop cluster kicked off with whirr.

Usage:
    cluster_install_distblast.py <cluster config file> <private_key_file>
"""
import os
import sys
import subprocess

import fabric.api as fabric
import fabric.contrib.files as fabric_files

def main(cluster_config, key_file):
    if cluster_config.endswith(".properties"):
        addresses = _get_whirr_addresses(cluster_config)
    else:
        addresses = _get_python_addresses(cluster_config)
    for addr in addresses:
        install_distblast(addr, key_file)

def install_distblast(addr, key_file):
    print "Installing on", addr
    with fabric.settings(host_string="%s@%s" % ("ubuntu", addr),
                         key_filename=key_file):
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

def _get_python_addresses(cluster_name):
    """Retrieve machine addresses using the older python hadoop-ec2 scripts.
    """
    cl = ["hadoop-ec2", "list", cluster_name]
    return _addresses_from_cl(cl)

def _get_whirr_addresses(whirr_config):
    """Retrieve IP addresses of cluster machines from Whirr.
    """
    cl = ["whirr", "list-cluster", "--config", whirr_config]
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
