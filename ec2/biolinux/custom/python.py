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
                run("python setup.py build")
                sudo("python setup.py install --skip-build")
                sudo("rm -rf dist")
                sudo("rm -rf lib/bx_python.egg-info")
