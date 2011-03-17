"""Install instructions for python libraries not ready for easy_install.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir, _if_not_python_lib

@_if_not_python_lib("bx")
def install_bx_python(env):
    """Install bx-python
    """
    clone_url = "http://bitbucket.org/james_taylor/bx-python"
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("hg clone %s" % clone_url)
            with cd(os.path.split(clone_url)[-1]):
                run("python%s setup.py build" % env.python_version_ext)
                sudo("python%s setup.py install --skip-build" % env.python_version_ext)
                sudo("rm -rf dist")
                sudo("rm -rf lib/bx_python.egg-info")
