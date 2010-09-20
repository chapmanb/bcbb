"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir

def install_pydoop(env):
    """Install pydoop; provides Hadoop access for Python.

    http://pydoop.sourceforge.net/docs/installation.html
    """
    with settings(warn_only=True):
        result = run("python -c 'import pydoop'")
    if not result.failed:
        return
    hadoop_version = "0.20.2"
    pydoop_version = "0.3.6"
    hadoop_url = "http://apache.thelorne.com/hadoop/core/" \
            "hadoop-%s/hadoop-%s.tar.gz" % (hadoop_version, hadoop_version)
    pydoop_url ="http://downloads.sourceforge.net/project/pydoop/" \
            "Pydoop/pydoop-%s.tar.gz" % (pydoop_version)

    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % hadoop_url)
            run("tar -xzvpf %s" % os.path.split(hadoop_url)[-1])
            hadoop_dir = os.path.join(work_dir, "hadoop-%s" % hadoop_version)
            run("wget %s" % pydoop_url)
            run("tar -xzvpf %s" % os.path.split(pydoop_url)[-1])
            with cd("pydoop-%s" % pydoop_version):
                run("export HADOOP_HOME=%s && export JAVA_HOME=%s && "\
                    "python setup.py build" % (hadoop_dir, env.java_home))
                sudo("python setup.py install --skip-build")

def install_mahout(env):
    # ToDo setup mahout, must be checked out from repo ATM:
    # https://cwiki.apache.org/MAHOUT/mahoutec2.html
    pass
