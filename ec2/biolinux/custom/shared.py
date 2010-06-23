"""Reusable decorators and functions for custom installations.
"""
import os
from contextlib import contextmanager

from fabric.api import *
from fabric.contrib.files import *

def _if_not_installed(pname):
    """Decorator that checks if a callable program is installed.
    """
    def argcatcher(func):
        def decorator(*args, **kwargs):
            with settings(
                    hide('warnings', 'running', 'stdout', 'stderr'),
                    warn_only=True):
                result = run(pname)
            if result.return_code == 127:
                return func(*args, **kwargs)
        return decorator
    return argcatcher

@contextmanager
def _make_tmp_dir():
    home_dir = run("echo $HOME")
    work_dir = os.path.join(home_dir, "tmp")
    if not exists(work_dir):
        run("mkdir %s" % work_dir)
    yield work_dir
    if exists(work_dir):
        run("rm -rf %s" % work_dir)
