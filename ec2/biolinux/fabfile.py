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
import sys
import subprocess

from fabric.api import *
from fabric.contrib.files import *
import yaml

env.config_dir = os.path.join(os.path.dirname(__file__), "config")

def ec2_ubuntu_environment():
    """Setup default environmental variables for Ubuntu EC2 servers.

    Works on a US EC2 server running Ubuntu. This is fairly general but
    we will need to define similar functions for other targets.
    """
    if env.hosts == ["vagrant"]:
        _setup_vagrant_environment()
    else:
        env.user = "ubuntu"
    env.sources_file = "/etc/apt/sources.list"
    env.shell_config = "~/.bashrc"
    # Global installation directory for packages and standard programs
    env.system_install = "/usr"
    # Local install directory for versioned software that will not
    # be included in the path by default.
    env.local_install = "install"
    env.shell = "/bin/bash -l -c"
    # XXX look for a way to find JAVA_HOME automatically
    env.java_home = "/usr/lib/jvm/java-6-openjdk"
    version = ("lucid", "10.4")
    #version = ("maverick", "10.10")
    sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb http://archive.canonical.com/ %s partner",
      "deb http://downloads.mongodb.org/distros/ubuntu % s 10gen",
      "deb http://watson.nci.nih.gov/cran_mirror/bin/linux/ubuntu %s/",
      "deb http://cran.stat.ucla.edu/bin/linux/ubuntu %s/", # lastest R versions
      "deb http://nebc.nox.ac.uk/bio-linux/ unstable bio-linux", # Bio-Linux
      "deb http://archive.cloudera.com/debian %s-cdh3 contrib", # Hadoop
      "ppa:freenx-team/ppa", # FreeNX PPA
      "deb http://download.virtualbox.org/virtualbox/debian %s contrib", # virtualbox
    ]
    env.std_sources = _add_source_versions(version, sources)

def _setup_vagrant_environment():
    """Use vagrant commands to get connection information.
    https://gist.github.com/1d4f7c3e98efdf860b7e
    """
    raw_ssh_config = subprocess.Popen(["vagrant", "ssh-config"],
                                      stdout=subprocess.PIPE).communicate()[0]
    ssh_config = dict([l.strip().split() for l in raw_ssh_config.split("\n") if l])
    env.user = ssh_config["User"]
    env.hosts = [ssh_config["HostName"]]
    env.port = ssh_config["Port"]
    env.host_string = "%s@%s:%s" % (env.user, env.hosts[0], env.port)
    env.key_filename = ssh_config["IdentityFile"]

def _add_source_versions(version, sources):
    name, num = version
    final = []
    for s in sources:
        if s.find("%s") > 0:
            s = s % name
        elif s.find("% s") > 0:
            s = s % num
        final.append(s)
    return final

def install_biolinux():
    """Main entry point for installing Biolinux on a remote server.
    """
    ec2_ubuntu_environment()
    pkg_install, lib_install = _read_main_config()
    _setup_sources()
    _setup_automation()
    _add_gpg_keys()
    _apt_packages(pkg_install)
    _custom_installs(pkg_install)
    _do_library_installs(lib_install)
    _freenx_scripts()
    _cleanup()

def _apt_packages(to_install):
    """Install packages available via apt-get.
    """
    pkg_config = os.path.join(env.config_dir, "packages.yaml")
    sudo("apt-get update")
    sudo("apt-get -y --force-yes upgrade")
    # Retrieve packages to get and install each of them
    (packages, _) = _yaml_to_packages(pkg_config, to_install)
    for package in packages:
        sudo("apt-get -y --force-yes install %s" % package)

def _custom_installs(to_install):
    if not exists(env.local_install):
        run("mkdir %s" % env.local_install)
    pkg_config = os.path.join(env.config_dir, "custom.yaml")
    packages, pkg_to_group = _yaml_to_packages(pkg_config, to_install)
    sys.path.append(os.path.split(__file__)[0])
    for p in packages:
        install_custom(p, True)

def install_custom(p, automated=False):
    """Install a single custom package by name.

    fab install_custom_package:package_name
    """
    if not automated:
        pkg_config = os.path.join(env.config_dir, "custom.yaml")
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)
        sys.path.append(os.path.split(__file__)[0])
        env.system_install = "/usr"
    try:
        mod = __import__("custom.%s" % pkg_to_group[p], fromlist=["custom"])
    except ImportError:
        raise ImportError("Need to write a %s module in custom." %
                pkg_to_group[p])
    try:
        fn = getattr(mod, "install_%s" % p)
    except AttributeError:
        raise ImportError("Need to write a install_%s function in custom.%s"
                % (p, pkg_to_group[p]))
    fn(env)

def _yaml_to_packages(yaml_file, to_install):
    """Read a list of packages from a nested YAML configuration file.
    """
    # allow us to check for packages only available on 64bit machines
    machine = run("uname -m")
    is_64bit = machine.find("_64") > 0
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    # filter the data based on what we have configured to install
    data = [(k, v) for (k,v) in full_data.iteritems()
            if to_install is None or k in to_install]
    data.sort()
    packages = []
    pkg_to_group = dict()
    while len(data) > 0:
        cur_key, cur_info = data.pop(0)
        if cur_info:
            if isinstance(cur_info, (list, tuple)):
                packages.extend(sorted(cur_info))
                for p in cur_info:
                    pkg_to_group[p] = cur_key
            elif isinstance(cur_info, dict):
                for key, val in cur_info.iteritems():
                    # if we are okay, propagate with the top level key
                    if key != 'needs_64bit' or is_64bit:
                        data.append((cur_key, val))
            else:
                raise ValueError(cur_info)
    return packages, pkg_to_group

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
    return packages, sorted(libraries)

# --- Library specific installation code

def _r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.
    """
    # Create an Rscript file with install details.
    out_file = "install_packages.R"
    run("touch %s" % out_file)
    repo_info = """
    cran.repos <- getOption("repos")
    cran.repos["CRAN" ] <- "%s"
    options(repos=cran.repos)
    source("%s")
    """ % (config["cranrepo"], config["biocrepo"])
    append(repo_info, out_file)
    install_fn = """
    repo.installer <- function(repos, install.fn) {
      update.or.install <- function(pname) {
        if (pname %in% installed.packages())
          update.packages(lib.loc=c(pname), repos=repos, ask=FALSE)
        else
          install.fn(pname)
      }
    }
    """
    append(install_fn, out_file)
    std_install = """
    std.pkgs <- c(%s)
    std.installer = repo.installer(cran.repos, install.packages)
    lapply(std.pkgs, std.installer)
    """ % (", ".join('"%s"' % p for p in config['cran']))
    append(std_install, out_file)
    bioc_install = """
    bioc.pkgs <- c(%s)
    bioc.installer = repo.installer(biocinstallRepos(), biocLite)
    lapply(bioc.pkgs, bioc.installer)
    """ % (", ".join('"%s"' % p for p in config['bioc']))
    append(bioc_install, out_file)
    final_update = """
    update.packages(repos=biocinstallRepos(), ask=FALSE)
    update.packages(ask=FALSE)
    """
    append(final_update, out_file)
    # run the script and then get rid of it
    sudo("Rscript %s" % out_file)
    run("rm -f %s" % out_file)

def _python_library_installer(config):
    """Install python specific libraries using easy_install.
    """
    for pname in config['pypi']:
        sudo("easy_install -U %s" % pname)

def _ruby_library_installer(config):
    """Install ruby specific gems.
    """
    for gem in config['gems']:
        with settings(
                hide('warnings', 'running', 'stdout', 'stderr')):
            gem_info = run("gem list --no-versions")
        installed = [l for l in gem_info.split("\n") if l]
        if gem in installed:
            sudo("gem update %s" % gem)
        else:
            sudo("gem install %s" % gem)

def _perl_library_installer(config):
    """Install perl libraries from CPAN with cpanminus.
    """
    run("wget --no-check-certificate http://xrl.us/cpanm")
    run("chmod a+rwx cpanm")
    sudo("mv cpanm %s/bin" % env.system_install)
    for lib in config['cpan']:
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        sudo("cpanm --skip-installed %s < /dev/null" % (lib))

def _clojure_library_installer(config):
    """Install clojure libraries using cljr.
    """
    for lib in config['cljr']:
        run("cljr install %s" % lib)

lib_installers = {
        "r-libs" : _r_library_installer,
        "python-libs" : _python_library_installer,
        "ruby-libs" : _ruby_library_installer,
        "perl-libs" : _perl_library_installer,
        "clojure-libs": _clojure_library_installer,
        }

def _do_library_installs(to_install):
    for iname in to_install:
        yaml_file = os.path.join(env.config_dir, "%s.yaml" % iname)
        with open(yaml_file) as in_handle:
            config = yaml.load(in_handle)
        lib_installers[iname](config)

# --- System hacks to support automation on apt systems

def _add_gpg_keys():
    """Adds GPG keys from all repositories
    """
    standalone = [
        "http://archive.cloudera.com/debian/archive.key",
        "http://download.virtualbox.org/virtualbox/debian/oracle_vbox.asc"]
    keyserver = [
        ("subkeys.pgp.net", "7F0CEB10"),
        ("subkeys.pgp.net", "381BA480"),
        ("keyserver.ubuntu.com", "D67FC6EAE2A11821"),
    ]
    for url, key in keyserver:
        sudo("apt-key adv --keyserver %s --recv %s" % (url, key))
    for key in standalone:
        sudo("wget -q -O- %s | apt-key add -" % key)

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
            "grub-pc grub2/linux_cmdline string ''",
            "grub-pc grub-pc/install_devices_empty boolean true",
            "acroread acroread/default-viewer boolean false",
            ]
    for l in package_info:
        sudo("echo %s | /usr/bin/debconf-set-selections" % l)

def _setup_sources():
    """Add sources for retrieving library packages.
       Using add-apt-repository allows processing PPAs
    """
    sudo("apt-get install -y --force-yes python-software-properties")
    for source in env.std_sources:
        if source.startswith("ppa:"):
            sudo("add-apt-repository '%s'" % source)
        elif not contains(source, env.sources_file):
            append(source, env.sources_file, use_sudo=True)

def _freenx_scripts():
    """Provide graphical access to clients via FreeNX.
    """
    setup_script = "setupnx.sh"
    remote_setup = "%s/bin/%s" % (env.system_install, setup_script)
    install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
    if not exists(remote_setup):
        put(os.path.join(install_file_dir, setup_script), setup_script,
                mode=0777)
        sudo("mv %s %s" % (setup_script, remote_setup))
    remote_login = "~/.bash_login"
    if not exists(remote_login):
        put(os.path.join(install_file_dir, 'bash_login'), remote_login,
                mode=0777)

def _cleanup():
    """Clean up any extra files after building.
    """
    run("rm -f $HOME/.bash_history")
    sudo("rm -f /var/crash/*")
    sudo("rm -f /var/log/firstboot.done")
    sudo("rm -f $HOME/.nx_setup_done")
    # RabbitMQ fails to start if its database is embedded into the image
    # because it saves the current IP address or host name so delete it now.
    # When starting up, RabbitMQ will recreate that directory.
    sudo('/etc/init.d/rabbitmq-server stop')
    for db_location in ['/var/lib/rabbitmq/mnesia', '/mnesia']:
        if exists(db_location):
            sudo('rm -rf %s' % db_location)
