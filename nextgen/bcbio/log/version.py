"""
Functions for getting version strings from software

The default approach for getting the version string from a program is to execute the 
program with the parameter '-v' and parse the output for a 'Version: XXXXX' string.
If this is not the correct approach for a particular piece of software, a custom 
regular expression that will capture the version string and/or a custom set of parameters 
that will be passed to the program in order to get it to print out the version string can 
be specified. If this approach cannot be used, custom routines can be defined and are
expected to return the version string. In addition, if the program is a script that is 
part of the bcbb pipeline, the current commit hash will be read from a file in the user's
home directory that was written during install.
"""

import os
import subprocess
import shlex
import re
import glob
import pkg_resources
from sys import stderr

# Regexp's defining how to extract the software version from each program's output
# The keys in the dictionary should match the keys under the 'program' section in
# the post-process.yaml configuration file. The regexp should contain one capturing 
# group that will capture the version string.
generic_regexp = r'Version:\s*(\S+)'
regexp = {
          'pdflatex': r'(3\.14.*)',
          'bowtie': r'bowtie\s+version\s+(\S+)',
          'fastqc': r'v(\S+)',
          'gatk': r'\(GATK\)\s+v(\S+?),',
          'ucsc_bigwig': r'wigToBigWig\s+v\s+(.+?)\s+\-',
          'picard': r'picard\-(\S+)'
          }

# The parameters to pass to each program in order for it to print the version. By default
# the program will be called with a -v option so if it needs to be invoked without arguments,
# specify an empty string here.
generic_parameter = '-v'
parameter = {
             'maq': '',
             'samtools': '',
             'bowtie': '--version',
             'bwa': '',
             'ucsc_bigwig': ''
             }

# The name of a file in the user's home directory where the version of the bcbb pipeline was written in setup.py
pipeline_version_file = '.bcbb_pipeline_version'
    
# A dictionary pointing to routines that know how to get the version info from a piece of software and overrides the generic version fetch
templates = {
             'picard': '_picard_version',
             'gatk': '_gatk_version'
             }

def export_git_commit_hash():
    v = _generic_version("git","rev-parse --short --verify HEAD",r'([a-d0-9]+)')
    if v != 'N/A':
        try:
            versionfile = os.path.join(os.getenv('HOME'),pipeline_version_file)
            with open(versionfile,'w') as vf:
                vf.write(v)
        except Exception, e:
            stderr.write("%s" % e)
            return 1
    return 0
    
def _generic_version(exe,param,regx):
    
    """Extract the version number of a piece of software.
    
    argument -- a string containing the arguments needed for the software to print the version number
    regexp -- a regular expression that will capture the version string from the output produced. The regexp must contain at least one capturing group and the first group is assumed to contain the version string
    
    """
    if param is None:
        param = generic_parameter
    if regx is None:
        regx = generic_regexp
        
    try:
        output = _get_output(exe,param)
    except Exception, e:
        log.info("Exception caught when trying to determine software version for %s: %s" % (exe,e))
        output = ""
        
    m = re.search(regx,output,re.I)
    if m and len(m.groups()) > 0:
        return m.group(1)
    return "N/A"

def _get_git_commit():
    
    v = 'N/A'
    # Look for a pipeline version file in the user's home directory
    versionfile = os.path.join(os.getenv('HOME'),pipeline_version_file)    
    if os.path.exists(versionfile):
        with open(versionfile) as vf:
            v = vf.read().strip()
    return v 

def _get_output(exe,param=''):
    args = shlex.split("%s %s" % (exe,param))
    output = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
    return output
    
def _get_pipeline_version():
    v = 'N/A'
    try:
        dist = pkg_resources.get_distribution("bcbio-nextgen")
    except:
        pass
    if dist:
        v = dist.version
    return v
 
def get_version(name,exe):
 
    """Get the software version of a named program and the corresponding executable"""
    
    name = name.lower()
    param = parameter.get(name,None)
    regx = regexp.get(name,None)
    
    # If the executable is a python script, attempt to get the path and the git commit hash for the repository
    if os.path.splitext(exe)[1] == '.py':
        version = _get_pipeline_version()
        
    # Else, if we have a function that knows how to get the version info from the particular piece of software, call that
    elif name in templates:
        version = globals()[templates[name]](name,exe,param,regx)
        
    # Else, try the generic approach..
    else:
        version = _generic_version(exe,param,regx)
        
    return version
    
def get_versions(config):
    
    """Returns a dictionary with the program names (as specified in the 'program' section of the supplied configuration object) as keys and the corresponding software version strings as values. If no version could be determined, the string 'N/A' is used."""
       
    # Get the list of external software from the config file
    prog_version = dict()
    for name, executable in config.get("program",{}).items():
        
        prog_version[name] = get_version(name,executable)
        
    return prog_version
       
def _gatk_version(name,path,param,regx):
    return _generic_version('java',"-jar %s --help" % os.path.join(path,"GenomeAnalysisTK.jar"),regx)    
def _picard_version(name,path,param,regx):
    # Assume that the picard jarfile has the version number in the file name...
    jars = glob.glob(os.path.join(path,'picard-*.jar'))
    if len(jars) > 0:    
        jarname = os.path.splitext(os.path.split(jars[0])[1])[0]
        m = re.search(r'picard\-(\S+)',jarname,re.I)
        if m and len(m.groups()) > 0:
            return m.group(1)
    return "N/A"
