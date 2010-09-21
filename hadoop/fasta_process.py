#!/usr/bin/env python
"""Process a fasta file through Hadoop one record at a time.
"""
import sys, os, logging, struct
import tempfile
import subprocess
import contextlib
logging.basicConfig(level=logging.DEBUG)

from pydoop.pipes import Mapper, Reducer, Factory, runTask
from pydoop.pipes import RecordReader, InputSplit, RecordWriter
from pydoop.hdfs import hdfs
from pydoop.utils import split_hdfs_path

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def top_blast_result(fasta_str, config):
    """Retrieve the top BLAST hit for a fasta record.
    """
    tmp_dir = config.get("job.local.dir")
    xref_db = config.get("fasta.blastdb")
    with tmpfile(prefix="in", dir=tmp_dir) as input_ref:
        with tmpfile(prefix="out", dir=tmp_dir) as blast_out:
            with open(input_ref, "w") as out_handle:
                out_handle.write(fasta_str)

            cl = NcbiblastnCommandline(query=input_ref, db=xref_db,
                    out=blast_out, outfmt=5, num_descriptions=1,
                    num_alignments=0)
            subprocess.check_call(str(cl).split())
            with open(blast_out) as blast_handle:
                rec = NCBIXML.read(blast_handle)
            xref = (rec.descriptions[0].title if len(rec.descriptions) > 0
                    else "nohit")
    return str(xref)

class FastaMapper(Mapper):
    def map(self, context):
        jc = context.getJobConf()
        result = top_blast_result(context.getInputValue(), jc)
        context.emit(context.getInputKey(), result)

class FastaReducer(Reducer):
    """Simple reducer that returns a value per input record identifier.
    """
    def reduce(self, context):
        key = context.getInputKey()
        vals = []
        while context.nextValue():
            vals.append(context.getInputValue())
        if len(vals) > 0:
            context.emit(key, vals[0])

class FastaReader(RecordReader):
    """Return one text FASTA record at a time using Biopython SeqIO iterators.
    """
    def __init__(self, context):
        super(FastaReader, self).__init__()
        self.logger = logging.getLogger(self.__class__.__name__)
        self.isplit = InputSplit(context.getInputSplit())
        self.host, self.port, self.fpath = split_hdfs_path(self.isplit.filename)
        self.fs = hdfs(self.host, self.port)
        self.file = self.fs.open_file(self.fpath, os.O_RDONLY)
        self._iterator = (SeqIO.parse(self.file, "fasta") if
                          self.isplit.offset == 0 else None)

    def __del__(self):
        self.file.close()
        self.fs.close()

    def next(self):
        if self._iterator:
            try:
                record = self._iterator.next()
                return (True, record.id, record.format("fasta"))
            except StopIteration:
                pass
        return (False, "", "")

    def getProgress(self):
        return 0

@contextlib.contextmanager
def tmpfile(*args, **kwargs):
    """Make a tempfile, safely cleaning up file descriptors on completion.
    """
    (fd, fname) = tempfile.mkstemp(*args, **kwargs)
    try:
        yield fname
    finally:
        os.close(fd)
        os.unlink(fname)

def main(argv):
    runTask(Factory(FastaMapper, FastaReducer,
                  record_reader_class=FastaReader,
                  ))

if __name__ == "__main__":
    main(sys.argv)
