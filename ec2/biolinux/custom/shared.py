"""Reusable decorators and functions for custom installations.
"""
import os
from contextlib import contextmanager

from fabric.api import *
from fabric.contrib.files import *

# -- decorators and context managers

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

def _if_not_python_lib(library):
    """Decorator that checks if a python library is installed.
    """
    def argcatcher(func):
        def decorator(*args, **kwargs):
            with settings(warn_only=True):
                result = run("python -c 'import %s'" % library)
            if result.failed:
                return func(*args, **kwargs)
        return decorator
    return argcatcher

@contextmanager
def _make_tmp_dir():
    home_dir = run("echo $HOME")
    work_dir = os.path.join(home_dir, "tmp", "cloudbiolinux")
    if not exists(work_dir):
        run("mkdir -p %s" % work_dir)
    yield work_dir
    if exists(work_dir):
        run("rm -rf %s" % work_dir)

# -- Standard build utility simplifiers

def _get_expected_file(url):
    tar_file = os.path.split(url)[-1]
    exts = {(".tar.gz", ".tgz") : "tar -xzpf",
            (".tar.bz2",): "tar -xjpf",
            (".zip",) : "unzip"}
    for ext_choices, tar_cmd in exts.iteritems():
        for ext in ext_choices:
            if tar_file.endswith(ext):
                return tar_file, tar_file[:-len(ext)], tar_cmd
    raise ValueError("Did not find extract command for %s" % url)

def _safe_dir_name(dir_name):
    replace_try = ["", "-src"]
    for replace in replace_try:
        check = dir_name.replace(replace, "")
        if exists(check):
            return check
    # still couldn't find it, it's a nasty one
    first_part = dir_name.split("-")[0].split("_")[0]
    with settings(warn_only=True):
        dirs = run("ls -d1 *%s*/" % first_part).split("\n")
    if len(dirs) == 1:
        return dirs[0]
    raise ValueError("Could not find directory %s" % dir_name)

def _fetch_and_unpack(url):
    if url.startswith(("git", "svn", "hg", "cvs")):
        run(url)
        base = os.path.basename(url.split()[-1])
        return os.path.splitext(base)[0]
    else:
        run("wget %s" % url)
        tar_file, dir_name, tar_cmd = _get_expected_file(url)
        run("%s %s" % (tar_cmd, tar_file))
        return _safe_dir_name(dir_name)

def _configure_make(env):
    run("./configure --prefix=%s " % env.system_install)
    run("make")
    sudo("make install")

def _make_copy(find_cmd, premake_cmd=None, do_make=True):
    def _do_work(env):
        install_dir = os.path.join(env.system_install, "bin")
        if premake_cmd:
            premake_cmd()
        if do_make:
            run("make")
        for fname in run(find_cmd).split("\n"):
            sudo("mv -f %s %s" % (fname, install_dir))
    return _do_work

def _get_install(url, env, make_commands):
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            dir_name = _fetch_and_unpack(url)
            with cd(dir_name):
                make_commands(env)

# --- Language specific utilities

def _symlinked_java_version_dir(pname, version):
    base_dir = os.path.join(env.system_install, "share", "java", pname)
    install_dir = "%s-%s" % (base_dir, version)
    if not exists(install_dir):
        sudo("mkdir -p %s" % install_dir)
        if exists(base_dir):
            sudo("rm -f %s" % base_dir)
        sudo("ln -s %s %s" % (install_dir, base_dir))
        return install_dir
    return None
