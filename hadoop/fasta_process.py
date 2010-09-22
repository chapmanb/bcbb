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


from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def top_blast_result(cur_id, fasta_str, xref_dbs, tmp_dir):
    """Retrieve the top BLAST hit for a fasta record.
    """
    xrefs = [cur_id]
    scores = [cur_id]
    with tmpfile(prefix="in", dir=tmp_dir) as input_ref:
        for xref_db in xref_dbs:
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
                xrefs.append(str(xref))
    return xrefs, scores

class FastaMapper(Mapper):
    def map(self, context):
        jc = context.getJobConf()
        tmp_dir = config.get("job.local.dir")
        xref_dbs = config.get("fasta.blastdb").split(",")
        ids, scores = blast.blast_top_hits(context.getInputKey(),
                context.getInputValue(), xref_dbs, tmp_dir)
        cur_key = ids[0]
        cur_val = dict(ids=ids[1:], scores=scores[1:])
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
