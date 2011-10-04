"""Helpful utilities for building analysis pipelines.
"""
import os
import tempfile
import shutil
import contextlib
import itertools
import functools
import ConfigParser
import csv, codecs, cStringIO
import subprocess
import shlex
import re
import glob
from sys import stderr

try:
    import multiprocessing
    from multiprocessing.pool import IMapIterator
except ImportError:
    multiprocessing = None

@contextlib.contextmanager
def cpmap(cores=1):
    """Configurable parallel map context manager.

    Returns appropriate map compatible function based on configuration:
    - Local single core (the default)
    - Multiple local cores
    """
    if int(cores) == 1:
        yield itertools.imap
    else:
        if multiprocessing is None:
            raise ImportError("multiprocessing not available")
        # Fix to allow keyboard interrupts in multiprocessing: https://gist.github.com/626518
        def wrapper(func):
            def wrap(self, timeout=None):
                return func(self, timeout=timeout if timeout is not None else 1e100)
            return wrap
        IMapIterator.next = wrapper(IMapIterator.next)
        pool = multiprocessing.Pool(int(cores))
        yield pool.imap
        pool.terminate()

def map_wrap(f):
    """Wrap standard function to easily pass into 'map' processing.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return apply(f, *args, **kwargs)
    return wrapper

def memoize_outfile(ext):
    """Creates outfile from input file and ext, running if outfile not present.

    This requires a standard function usage. The first arg, or kwarg 'in_file', needs
    to be the input file that is being processed. The output name is created with the
    provided ext relative to the input. The function is only run if the created
    out_file is not present.
    """
    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            if len(args) > 0:
                in_file = args[0]
            else:
                in_file = kwargs['in_file']
            out_file = "%s%s" % (os.path.splitext(in_file)[0], ext)
            if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                kwargs['out_file'] = out_file
                f(*args, **kwargs)
            return out_file
        return wrapper
    return decor

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if not os.path.isdir(dname):
                raise
    return dname

@contextlib.contextmanager
def curdir_tmpdir(remove=True):
    """Context manager to create and remove a temporary directory.
    """
    tmp_dir_base = os.path.join(os.getcwd(), "tmp")
    safe_makedir(tmp_dir_base)
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_base)
    safe_makedir(tmp_dir)
    try :
        yield tmp_dir
    finally :
        if remove:
            shutil.rmtree(tmp_dir)

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    safe_makedir(new_dir)
    os.chdir(new_dir)
    try :
        yield
    finally :
        os.chdir(cur_dir)

@contextlib.contextmanager
def tmpfile(*args, **kwargs):
    """Make a tempfile, safely cleaning up file descriptors on completion.
    """
    (fd, fname) = tempfile.mkstemp(*args, **kwargs)
    try:
        yield fname
    finally:
        os.close(fd)
        if os.path.exists(fname):
            os.remove(fname)

@contextlib.contextmanager
def file_transaction(*rollback_files):
    """Wrap file generation in a transaction, removing partial files on failure.

   This allows a safe restart at any point, helping to deal with interrupted
   pipelines.
   """
    try:
        yield None
    except:
        for fnames in rollback_files:
            if isinstance(fnames, str):
                fnames = [fnames]
            for fname in fnames:
                if fname and os.path.exists(fname) and os.path.isfile(fname):
                    os.remove(fname)
        raise

def create_dirs(config, names=None):
    if names is None:
        names = config["dir"].keys()
    for dname in names:
        d = config["dir"][dname]
        safe_makedir(d)

def save_diskspace(fname, reason, config):
    """Overwrite a file in place with a short message to save disk.

    This keeps files as a sanity check on processes working, but saves
    disk by replacing them with a short message.
    """
    if config["algorithm"].get("save_diskspace", False):
        with open(fname, "w") as out_handle:
            out_handle.write("File removed to save disk space: %s" % reason)

def read_galaxy_amqp_config(galaxy_config, base_dir):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    galaxy_config = add_full_path(galaxy_config, base_dir)
    config = ConfigParser.ConfigParser()
    config.read(galaxy_config)
    amqp_config = {}
    for option in config.options("galaxy_amqp"):
        amqp_config[option] = config.get("galaxy_amqp", option)
    return amqp_config

def add_full_path(dirname, basedir=None):
    if basedir is None:
        basedir = os.getcwd()
    if not dirname.startswith("/"):
        dirname = os.path.join(basedir, dirname)
    return dirname

# UTF-8 methods for csv module (does not support it in python >2.7)
# http://docs.python.org/library/csv.html#examples

class UTF8Recoder:
    """Iterator that reads an encoded stream and reencodes the input to UTF-8
    """
    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def __iter__(self):
        return self

    def next(self):
        return self.reader.next().encode("utf-8")


class UnicodeReader:
    """A CSV reader which will iterate over lines in the CSV file "f",
       which is encoded in the given encoding.
    """
    
    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        f = UTF8Recoder(f, encoding)
        self.reader = csv.reader(f, dialect=dialect, **kwds)

    def next(self):
        row = self.reader.next()
        return [unicode(s, "utf-8") for s in row]

    def __iter__(self):
        return self

class UnicodeWriter:
    """A CSV writer which will write rows to CSV file "f",
       which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        # Redirect output to a queue
        self.queue = cStringIO.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def writerow(self, row):
        self.writer.writerow([str(s).encode("utf-8") for s in row])
        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()
        data = data.decode("utf-8")
        # ... and reencode it into the target encoding
        data = self.encoder.encode(data)
        # write to the target stream
        self.stream.write(data)
        # empty queue
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)

class SoftwareVersion:
    """A class that knows how to query external software used by the pipeline for
    the version information."""
    
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
    
    # A static method for writing the git commit hash to the pipeline_version_file
    @classmethod
    def export_git_commit_hash(cls):
        args = shlex.split("git rev-parse --short --verify HEAD")
        try:
            v = (subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].strip() or "N/A")
            # Only export the returned value if it is a hexidecimal string
            if re.search(r'[a-d0-9]+',v):
                versionfile = os.path.join(os.getenv('HOME'),cls.pipeline_version_file)
                with open(versionfile,'w') as vf:
                    vf.write(v)
        except Exception, e:
            stderr.write("%s" % e)
            return 1
        return 0
    
    def __init__(self,config):
        
        self.version = {}
        self.executable = {}
        
        # Parse the paths and executables in the supplied config and store them locally
        for name, exe in config.get("program",{}).items():
            name = name.lower()
            self.executable[name] = exe
            
            # For GATK, store the non-static data in the regexp and parameter dictionaries
            if name == 'gatk':
                self.executable[name] = 'java'
                self.parameter[name] = "-jar %s --help" % os.path.join(executable,"GenomeAnalysisTK.jar")
            
    def _get_executable(self,name):
        return self.executable.get(name,None)
    
    def _get_generic_version(self,name):
        
        v = None
        exe = self._get_executable(name)
        if exe:
            try:
                args = shlex.split("%s %s" % (exe,self._get_parameter(name)))
                output = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode('UTF-8')
                m = re.search(self._get_regexp(name),output,re.I)
                if m and len(m.groups()) > 0:
                    v = m.group(1)
            except Exception, e:
                stderr.write("Exception caught when trying to determine software version for %s: %s" % (exe,e))
        
        return v
    
    def _get_parameter(self,name):
        return self.parameter.get(name,self.generic_parameter)
    
    def _get_picard_version(self,name):
        
        v = None
        picard_home_dir = self._get_executable(name)
        if picard_home_dir:
            # Assume that the picard jarfile has the version number in the file name...
            jars = glob.glob(os.path.join(picard_home_dir,'picard-*.jar'))
            if len(jars) > 0:    
                jarname = os.path.splitext(os.path.split(jars[0])[1])[0]
                m = re.search(self._get_regexp(name),jarname,re.I)
                if m and len(m.groups()) > 0:
                    v = m.group(1)
                    
        return v

    def _get_pipeline_version(self):
        
        v = 'N/A'
        # Look for a pipeline version file in the user's home directory
        versionfile = os.path.join(os.getenv('HOME'),self.pipeline_version_file)    
        if os.path.exists(versionfile):
            with open(versionfile) as vf:
                v = vf.read().strip()
        return v 
   
    def _get_regexp(self,name):
        return self.regexp.get(name,self.generic_regexp)
            
    def version(self,name):
        """Get the version of the named software"""
        
        # Check if the version has been looked up and cached 
        name = name.lower()
        v = self.version(name,None)
        if v:
            return v
        
        # Check the special cases that cannot be queried generically:
        exe = self._get_executable(name)
        
        # Picard
        if name == 'picard':
            v = self._get_picard_version(name)
        # If the executable is a python script, get the pipeline version
        elif exe and os.path.splitext(exe)[1] == '.py':
            v = self._get_pipeline_version()
        else:
            v = self._get_generic_version(name)
        
        # cache the result and return
        self.version[name] = (v or 'N/A')
        return self.version[name]
    