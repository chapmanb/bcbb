"""
Misc paver tasks.

Several tasks rely on the options dictionary visible to all of pavers submodules. Make sure you have included the proper variables.
"""
import os

from paver.easy import *

@task
@cmdopts([('pattern=', 'p', 'pattern for rsync')])
def sync_to_dest():
    """Sync files in options.rsync.path with a pattern 'pat' to options.rsync.user@options.rsync.host:options.rsync.dest"""
    rsync = options.rsync
    opts = rsync.opts if rsync.has_key("opts") else "-av"
    path = os.path.realpath(rsync.path) + os.sep if rsync.has_key("path") else os.path.realpath(".") + os.sep
    dest = os.path.normpath(rsync.dest) + os.sep if rsync.has_key("dest") else os.path.normpath(".") + os.sep
    pat = options.sync_to_dest.pattern if options.sync_to_dest.has_key("pattern") else ""
    cl = [rsync.program, opts,  "--include='*/'", 
          "--include='" + pat + "'", "--exclude='*'",
          path, rsync.user + "@" + rsync.host + ":" + dest]
    sh(" ".join(cl))

@task
@cmdopts([('pattern=', 'p', 'pattern for rsync')])
def sync_from_dest():
    """Sync files in options.rsync.user@options.rsync.host:options.rsync.src with a pattern 'pat' to options.rsync.path"""
    rsync = options.rsync
    opts = rsync.opts if rsync.has_key("opts") else "-av"
    path = os.path.realpath(rsync.path) + os.sep if rsync.has_key("path") else os.path.realpath(".") + os.sep
    src = os.path.normpath(rsync.src) + os.sep if rsync.has_key("src") else os.path.normpath(".") + os.sep
    pat = options.sync_from_dest.pattern if options.sync_from_dest.has_key("pattern") else ""
    cl = [rsync.program, opts,  "--include='*/'", 
          "--include='" + pat + "'", "--exclude='*'",
          rsync.user + "@" + rsync.host + ":" + src, path]
    sh(" ".join(cl))

@task
def list_options():
    """List options"""
    print str(options)

### Functions for processing lists and dicts
def process_list(fn, items):
    """process list in function"""
    for item in items:
        fn(item)

def process_dict(fn, items):
    """process dict in function"""
    for s in items.keys():
        fn(items[s])

