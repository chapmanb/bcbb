"""Functions for getting version strings from software"""

import os
import subprocess
import shlex
import re
import glob

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
    
    # Look for a file named .[bcbb repository]_commit in the user's home directory
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
