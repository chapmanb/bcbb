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

env.config_dir = os.path.join(os.getcwd(), "config")

def ec2_ubuntu_environment():
    """Setup default environmental variables for Ubuntu EC2 servers.

    Works on a US EC2 server running Ubunutu 10.04 lucid. This should be pretty
    general but should support system specific things for other platform
    targets.
    """
    env.user = "ubuntu"
    env.sources_file = "/etc/apt/sources.list"
    env.std_sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ lucid universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ lucid universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ lucid-updates universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ lucid-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ lucid multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ lucid multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ lucid-updates multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ lucid-updates multiverse",
      "deb http://archive.canonical.com/ lucid partner",
    ]

def install_biolinux():
    """Main entry point for installing Biolinux on a remote server.
    """
    ec2_ubuntu_environment()
    pkg_install = _read_main_config()
    _apt_packages(pkg_install)

def _apt_packages(to_install):
    """Install packages available via apt-get.
    """
    pkg_config = os.path.join(env.config_dir, "packages.yaml")
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
    packages = _yaml_to_packages(pkg_config, to_install)
    _setup_automation()
    for package in packages:
        sudo("apt-get -y --force-yes install %s" % package)

def _setup_automation():
    """Setup the environment to be fully automated for installs.

    Sun Java license acceptance:
    http://www.davidpashley.com/blog/debian/java-license

    MySQL root password questions:
    http://snowulf.com/archives/540-Truly-non-interactive-unattended-apt-get-install.html
    """
    run("export DEBIAN_FRONTEND=noninteractive")
    license_info = [
            "sun-java6-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-bin shared/accepted-sun-dlj-v1-1 select true",
            ]
    for l in license_info:
        sudo("echo %s | /usr/bin/debconf-set-selections" % l)

def _yaml_to_packages(yaml_file, to_install):
    """Read a list of packages from a nested YAML configuration file.
    """
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    # filter the data based on what we have configured to install
    data = [(k, v) for (k,v) in full_data.iteritems() if k in to_install]
    data.sort()
    data = [v for (_, v) in data]
    packages = []
    while len(data) > 0:
        cur_info = data.pop(0)
        if cur_info:
            if isinstance(cur_info, (list, tuple)):
                packages.extend(sorted(cur_info))
            elif isinstance(cur_info, dict):
                data.extend(cur_info.values())
            else:
                raise ValueError(cur_info)
    return packages

def _read_main_config():
    """Pull a list of groups to install based on our main configuration YAML.
    """
    yaml_file = os.path.join(env.config_dir, "main.yaml")
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    return full_data['packages']
