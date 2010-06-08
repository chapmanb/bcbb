"""Main Fabric deployment file for BioLinux distribution.

This installs a standard set of useful biological applications on a remote
server. It is designed for bootstrapping a machine from scratch, as with new
Amazon EC2 instances.

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
    env.shell_config = "~/.bashrc"
    env.shell = "/bin/bash -l -c"
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
      "deb http://downloads.mongodb.org/distros/ubuntu 10.4 10gen",
      "deb http://archive.cloudera.com/debian karmic-cdh3b1 contrib",
      # ToDo packages for cloudera not available on lucid yet, using karmic for the moment (beta 1)
    ]

def install_biolinux():
    """Main entry point for installing Biolinux on a remote server.
    """
    ec2_ubuntu_environment()
    pkg_install, lib_install = _read_main_config()
    _add_gpg_keys()
    _apt_packages(pkg_install)
    _do_library_installs(lib_install)

def _add_gpg_keys():
    """Adds GPG keys from all repositories
       ToDo Cleanup/unify this
    """
    sudo("curl -s http://archive.cloudera.com/debian/archive.key | apt-key add -")
    # mongodb & CRAN
    sudo("apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10")
    sudo("gpg --keyserver subkeys.pgp.net --recv-key 381BA480")

def _apt_packages(to_install):
    """Install packages available via apt-get.
    """
    pkg_config = os.path.join(env.config_dir, "packages.yaml")
    # Setup and update apt sources on the remote host
    # lastest R versions and Bio-Linux. debian-med should already be there.
    sources_add = [
      "deb http://cran.stat.ucla.edu/bin/linux/ubuntu lucid/",
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
    """Setup the environment to be fully automated for tricky installs.

    Sun Java license acceptance:
    http://www.davidpashley.com/blog/debian/java-license

    MySQL root password questions; install with empty root password:
    http://snowulf.com/archives/540-Truly-non-interactive-unattended-apt-get-install.html

    Postfix, setup for no configuration. See more on issues here:
    http://www.uluga.ubuntuforums.org/showthread.php?p=9120196
    """
    interactive_cmd = "export DEBIAN_FRONTEND=noninteractive"
    if not contains(interactive_cmd, env.shell_config):
        append(interactive_cmd, env.shell_config)
    package_info = [
            "postfix postfix/main_mailer_type select No configuration",
            "postfix postfix/mailname string notusedexample.org",
            "mysql-server-5.1 mysql-server/root_password string '(password omitted)'",
            "mysql-server-5.1 mysql-server/root_password_again string '(password omitted)'",
            "sun-java6-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-bin shared/accepted-sun-dlj-v1-1 select true",
            ]
    for l in package_info:
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
    packages = full_data['packages']
    packages = packages if packages else []
    libraries = full_data['libraries']
    libraries = libraries if libraries else []
    return packages, libraries

# --- Library specific installation code

def _r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.

    ToDo: Can we check if a package is installed and only upgrade
    it if it's old? install.packages blindly downloads away.
    """
    # create an Rscript file which will be used to do the install
    out_file = "install_packages.R"
    run("touch %s" % out_file)
    # CRAN
    append(['repo <- getOption("repos")',
            'repo["CRAN" ] <- "%s"' % config["cranrepo"],
            'options(repos=repo)'], out_file)
    for pname in config['cran']:
        append('install.packages("%s")' % pname, out_file)
    # Bioconductor
    append('source("%s")' % config['biocrepo'], out_file)
    for pname in config['bioc']:
        append('biocLite("%s")' % pname, out_file)
    # run the script and then get rid of it
    sudo("Rscript %s" % out_file)
    run("rm -f %s" % out_file)

def _python_library_installer(config):
    """Install python specific libraries using easy_install.
    """
    for pname in config['pypi']:
        sudo("easy_install -U %s" % pname)

def _ruby_library_installer(config):
    """Install ruby specific gems using "gem install"
    """
    for gem in config['gems']:
	sudo("gem install %s" % gem)

# Note that the following Cloudera hadoop installation is for test
# purposes "inside the instance" only, Amazon already provides a
# production-ready Elastic MapReduce platform:
# http://aws.amazon.com/elasticmapreduce/

def _setup_hadoop(config):
    """Sets up Cloudera's CDH Hadoop and friends
	http://archive.cloudera.com/docs/ec2.html
	http://archive.cloudera.com/cdh/3/
    """

    # ToDo setup config files according to simple node config

    for pkg in config['mapreduce']:
    	sudo("apt-get install %s" % pkg)

    # ToDo setup mahout, must be checked out from repo ATM:
    # https://cwiki.apache.org/MAHOUT/mahoutec2.html

    #_checkout_repository()

    pass

#def _checkout_repository(url):
#     """ ToDo Checks out a repository
#     """
#    pass

lib_installers = {
        "r-libs" : _r_library_installer,
        "python-libs" : _python_library_installer,
	"ruby-libs" : _ruby_library_installer,
        }

def _do_library_installs(to_install):
    for iname in to_install:
        yaml_file = os.path.join(env.config_dir, "%s.yaml" % iname)
        with open(yaml_file) as in_handle:
            config = yaml.load(in_handle)
        lib_installers[iname](config)
