"""Script to extract the version of software that is run in the analysis pipeline.

This information is meant to be included in the report file and is crucial for
e.g. reproducing and publishing the results.

Usage:
    software_versions.py <post-processing config file>

The software used in the pipeline and which will be included in the report are specified in 
a post-processing configuration file, under the 'program' section. In case the program is a 
python script rather than an executable, the git commit hash will be reported.
"""
import os
import sys
import yaml
import logbook
import subprocess
import shlex
import re
from optparse import OptionParser

GIT_COMMIT_HASH=""

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

# A dictionary pointing to routines that know how to get the version info from a piece of software
templates = {
             'pdflatex': '_pdflatex_version'
             }

def main(config_file):

    # Parse the config yaml file
    if not os.path.exists(config_file):
        log.info("Could not find specified yaml configuration file %s." % config_file)
        return
    
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    
    # Get the list of currently loaded modules
    modules = dict()
    try:
        args = shlex.split("modulecmd bash -t list")
        output = subprocess.Popen(args, stderr=subprocess.PIPE).communicate()[1]
        rows = output.decode('UTF-8').split("\n")
        for row in output.decode('UTF-8').split("\n"):
            pcs = row.split("/")
            if (len(pcs) == 2):
                modules[pcs[0].lower()] = pcs[1]
    except Exception, e:
        log.error("Could not get module list: %s" % e)
    
    # Get the list of external software from the config file
    prog_version = dict()
    for name, executable in config.get("program",{}).items():
        
        name = name.lower()
        
        # If the executable is a python script, get the path and get the git commit hash for the repository
        if os.path.splitext(executable)[1] == '.py':
            version = _get_git_hash(executable)
        # Else, if the software is in the list of loaded modules, get the version from there
        elif name in modules:
            version = modules[name]
        # Else, if we have a function that knows how to get the version info from the particular piece of software, call that
        elif name in templates:
            version = globals()[templates[name]](executable)
        # Else, we don't know..
        else
            version = "N/A"
            
        prog_version[name] = version
        
    for name,version in prog_version.items():
        print "%s: version %s" % (name,version)
        
def _get_git_hash(script_name):
    
    if len(GIT_COMMIT_HASH) > 0:
        return GIT_COMMIT_HASH
    else:
        return "N/A"

def _generic_version(executable):
    
    args = shlex.split("%s -v" % executable)
    output = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
    version = output.decode('UTF-8')
    return version
    
def _pdflatex_version(executable):
    
    output = _generic_version(executable)
    m = re.search(r'(3\.14.*)',output.split("\n")[0])
    if m and len(m.groups()) > 0:
        return m.group(1)
    return "N/A"
         
if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)

