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
import glob

from optparse import OptionParser

# The name of the git repository (if any) where the python scripts are contained
git_repo = 'bcbb'

# A dictionary pointing to routines that know how to get the version info from a piece of software
templates = {
             'pdflatex': '_pdflatex_version',
             'maq': '_maq_version',
             'samtools': '_samtools_version',
             'bowtie': '_bowtie_version',
             'picard': '_picard_version',
             'bwa': '_bwa_version',
             'fastqc': '_fastqc_version',
             'gatk': '_gatk_version',
             'ucsc_bigwig': '_ucsc_bigwig_version'
             }

def main(config_file):

    # Parse the config yaml file
    if not os.path.exists(config_file):
        log.info("Could not find specified yaml configuration file %s." % config_file)
        return
    
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    for name,ver in get_versions(config).items():
        print "%s: %s" % (name,ver)

def get_versions(config):
    
    """Returns a dictionary with the program names (as specified in the 'program' section of the supplied configuration object) as keys and the corresponding software version strings as values. If no version could be determined, the string 'N/A' is used."""
       
    # Get the list of external software from the config file
    prog_version = dict()
    for name, executable in config.get("program",{}).items():
        
        prog_version[name] = get_version(name,executable)
        
    return prog_version

def get_version(name,executable):
 
    """Get the software version of a named program and the corresponding executable"""
    
    name = name.lower()
    
    # If the executable is a python script, attempt to get the path and the git commit hash for the repository
    if os.path.splitext(executable)[1] == '.py':
        version = _get_git_hash(git_repo)
        
    # Else, if we have a function that knows how to get the version info from the particular piece of software, call that
    elif name in templates:
        version = globals()[templates[name]](executable)
        
    # Else, try the generic approach..
    else:
        version = _generic_version(executable)
        
    return version
            
def _get_git_hash(git_repo):
    
    hash_file = os.path.join(os.getenv('HOME'),".%s_commit" % git_repo)
    commit = "N/A"
    if os.path.exists(hash_file):
        with open(hash_file) as hf:
            commit = hf.read().strip()
    return commit 
   
def _get_output(executable,argument=''):
    args = shlex.split("%s %s" % (executable,argument))
    output = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
    return output.decode('UTF-8')
     
def _generic_version(executable,argument='-v',regexp=r'Version:\s*(\S+)'):
    
    """Extract the version number of a piece of software.
    
    argument -- a string containing the arguments needed for the software to print the version number
    regexp -- a regular expression that will capture the version string from the output produced. The regexp must contain at least one capturing group and the first group is assumed to contain the version string
    
    """
    
    try:
        output = _get_output(executable,argument)
    except Exception, e:
        log.info("Exception caught when trying to determine software version for %s: %s" % (executable,e))
        output = ""
        
    m = re.search(regexp,output,re.I)
    if m and len(m.groups()) > 0:
        return m.group(1)
    return "N/A"

    
def _maq_version(executable):
    return _generic_version(executable,'')
def _samtools_version(executable):
    return _generic_version(executable,'')
def _bwa_version(executable):
    return _generic_version(executable,'')
def _ucsc_bigwig_version(executable):
    return _generic_version(executable,'',r'wigToBigWig\s+v\s+(.+?)\s+\-')
def _pdflatex_version(executable):
    return _generic_version(executable,'-v',r'(3\.14.*)') # The pdflatex version numbers are approaching pi...
def _fastqc_version(executable):
    return _generic_version(executable,'-v',r'v(\S+)')
def _bowtie_version(executable):
    return _generic_version(executable,'--version',r'bowtie\s+version\s+(\S+)')
def _gatk_version(executable):
    return _generic_version('java',"-jar %s --help" % os.path.join(executable,"GenomeAnalysisTK.jar"),r'\(GATK\)\s+v(\S+?),')    
def _picard_version(picard_home_dir):
    # Assume that the picard jarfile has the version number in the file name...
    jars = glob.glob(os.path.join(picard_home_dir,'picard-*.jar'))
    if len(jars) > 0:    
        jarname = os.path.splitext(os.path.split(jars[0])[1])[0]
        m = re.search(r'picard\-(\S+)',jarname,re.I)
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

