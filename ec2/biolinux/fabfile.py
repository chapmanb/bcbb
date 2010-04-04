"""Main Fabric deployment file for BioLinux distribution.

This installs a standard set of useful biological applications on a remote
server. It is designed for bootstrapping a machine from scratch, as with new
Amaxon EC2 instances.

Usage:
    fab -H hostname -i private_key_file install_biolinux

Requires:
    Fabric http://docs.fabfile.org
    PyYAML http://pyyaml.org/wiki/PyYAMLDocumentation
"""
import os

from fabric.api import *
from fabric.contrib.files import *
import yaml

def ec2_environment():
    """Setup default environmental variables for Ubuntu EC2 servers.
    """
    env.user = "root"
    env.sources_file = "/etc/apt/sources.list"

def install_biolinux():
    """Main entry point for installing Biolinux on a remote server.
    """
    ec2_environment()
    _apt_packages()

def _apt_packages():
    """Install packages available via apt-get.
    """
    pkg_config = os.path.join(os.getcwd(), "config", "packages.yaml")
    # Setup and update apt sources on the remote host
    # lastest R versions and Bio-Linux. debian-med should already be there.
    sources_add = [
      "deb http://cran.stat.ucla.edu/bin/linux/ubuntu karmic/",
      "deb http://nebc.nox.ac.uk/bio-linux/ unstable bio-linux",
      ]
    for source in sources_add:
        if not contains(source, env.sources_file):
            append(source, env.sources_file)
    run("apt-get update")
    # Retrieve packages to get and install each of them
    packages = _yaml_to_packages(pkg_config)
    for package in packages:
        run("apt-get -y install %s" % package)

def _yaml_to_packages(yaml_file):
    """Read a list of packages from a nested YAML configuration file.
    """
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    data = full_data.values()
    packages = []
    while len(data) > 0:
        cur_info = data.pop(0)
        if cur_info:
            if isinstance(cur_info, (list, tuple)):
                packages.extend(cur_info)
            elif isinstance(cur_info, dict):
                data.extend(cur_info.values())
            else:
                raise ValueError(cur_info)
    return packages
