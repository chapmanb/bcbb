#!/usr/bin/env python
"""Install distblast software on every node of a Hadoop cluster.

This is an example of how to remotely add non-AMI data or software
to a Hadoop cluster kicked off with whirr.

Usage:
    cluster_install_distblast.py <whirr config file>
"""
import os
import sys
import subprocess

import fabric.api as fabric
import fabric.contrib.files as fabric_files

def main(whirr_config):
    addresses = _get_cluster_addresses(whirr_config)
    for addr in addresses:
        install_distblast(addr)

def install_distblast(addr):
    print "Installing on", addr
    with fabric.settings(host_string="%s@%s" % ("ubuntu", addr)):
        work_dir = "install"
        if not fabric_files.exists(work_dir):
            fabric.run("mkdir %s" % work_dir)
        with fabric.cd(work_dir):
            fabric.run("git clone git://github.com/chapmanb/bcbb.git")
            with fabric.cd("bcbb/distblast"):
                fabric.run("python2.6 setup.py build")
                fabric.sudo("python2.6 setup.py install")

def _get_cluster_addresses(whirr_config):
    cl = ["whirr", "list-cluster", "--config", whirr_config]
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
