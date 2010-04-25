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

def ec2_ubuntu_environment():
    """Setup default environmental variables for Ubuntu EC2 servers.
    """
    env.user = "ubuntu"
    env.sources_file = "/etc/apt/sources.list"
    env.std_sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ karmic universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ karmic universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ karmic universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ karmic-updates universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ karmic-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ karmic multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ karmic multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ karmic-updates multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ karmic-updates multiverse",
    ]

def install_biolinux():
    """Main entry point for installing Biolinux on a remote server.
    """
    ec2_ubuntu_environment()
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
    for source in sources_add + env.std_sources:
        if not contains(source, env.sources_file):
            append(source, env.sources_file, use_sudo=True)
    sudo("apt-get update")
    # Retrieve packages to get and install each of them
    packages = _yaml_to_packages(pkg_config)
    _setup_licenses()
    for package in packages:
        sudo("apt-get -y --force-yes install %s" % package)

def _setup_licenses():
    """Handle automated license acceptance for things like Sun java.

    http://www.davidpashley.com/blog/debian/java-license
    """
    license_info = [
            "sun-java5-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java5-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-bin shared/accepted-sun-dlj-v1-1 select true",
            ]
    for l in license_info:
        sudo("echo %s | /usr/bin/debconf-set-selections" % l)

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
