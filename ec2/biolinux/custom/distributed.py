"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir, _if_not_python_lib, _fetch_and_unpack

@_if_not_python_lib("pydoop")
def install_pydoop(env):
    """Install pydoop; provides Hadoop access for Python.

    http://pydoop.sourceforge.net/docs/installation.html
    """
    hadoop_version = "0.21.0"
    pydoop_version = "0.3.7_rc1"
    hadoop_url = "http://apache.mirrors.hoobly.com/hadoop/core/" \
            "hadoop-%s/hadoop-%s.tar.gz" % (hadoop_version, hadoop_version)
    pydoop_url ="http://downloads.sourceforge.net/project/pydoop/" \
                "Pydoop-%s/pydoop-%s.tar.gz" % (pydoop_version, pydoop_version)

    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            hadoop_dir = _fetch_and_unpack(hadoop_url)
            pydoop_dir = _fetch_and_unpack(pydoop_url)
            with cd(pydoop_dir):
                export_str = "export HADOOP_HOME=%s && export JAVA_HOME=%s" % \
                    (os.path.join(os.pardir, hadoop_dir), env.java_home)
                run("%s && python%s setup.py build" % (export_str, env.python_version_ext))
                sudo("%s && python%s setup.py install --skip-build" %
                     (export_str, env.python_version_ext))

def install_mahout(env):
    # ToDo setup mahout, must be checked out from repo ATM:
    # https://cwiki.apache.org/MAHOUT/mahoutec2.html
    pass
