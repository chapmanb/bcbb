#!/usr/bin/env python
"""Process a fasta file through Hadoop one record at a time.
"""
import sys
import os
import json
import tempfile
import subprocess
import contextlib
import logging
logging.basicConfig(level=logging.DEBUG)

from pydoop.pipes import Mapper, Reducer, Factory, runTask
from pydoop.pipes import RecordReader, InputSplit, RecordWriter
from pydoop.hdfs import hdfs
from pydoop.utils import split_hdfs_path

from Bio import SeqIO

from bcbio.phylo import blast

class FastaMapper(Mapper):
    def map(self, context):
        config = context.getJobConf()
        tmp_dir = config.get("job.local.dir")
        xref_dbs = config.get("fasta.blastdb").split(",")
        cur_key, ids, scores = blast.blast_top_hits(context.getInputKey(),
                context.getInputValue(), xref_dbs, tmp_dir)
        cur_val = dict(ids=ids, scores=scores)
        context.emit(cur_key, json.dumps(cur_val))

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
