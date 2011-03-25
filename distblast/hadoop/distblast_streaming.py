#!/usr/bin/env python
"""Process a FASTA file using Hadoop streaming.

Examples with Dumbo and MrJob.
"""
import os
import sys
import json

from bcbio.phylo import blast

def mapper(key, rec):
    tmp_dir = os.environ["job_local_dir"]
    xref_dbs = os.environ["fasta_blastdb"].split(",")
    parts = rec.split("\t")
    if len(parts) == 3: # remove extra initial tab if present
        parts = parts[1:]
    title, seq = rec.split("\t")
    rec_id = title.split()[0]
    cur_key, ids, scores = blast.blast_top_hits(rec_id, seq, xref_dbs, tmp_dir)
    cur_val = dict(ids=ids, scores=scores)
    yield cur_key, cur_val

def reducer(key, vals):
    for v in vals:
        yield key, json.dumps(v)

# -- Alternative MrJob version.
want_mrjob=False

if want_mrjob:
    from mrjob.job import MRJob

    class DistblastJob(MRJob):
        def hadoop_job_runner_kwargs(self):
            config = MRJob.hadoop_job_runner_kwargs(self)
            # config["hadoop_extra_args"].extend(
            #         ["-inputformat", "com.lifetech.hadoop.streaming.FastaInputFormat",
            #          "-libjars", "jars/bioseq-0.0.1.jar"])
            return config

        def mapper(self, key, rec):
            for k, v in mapper(key, rec):
                yield k, v

        def reducer(self, key, vals):
            for k, v in reducer(key, vals):
                yield k, v

if __name__ == '__main__':
    if want_mrjob:
        DistblastJob.run()
    else:
        import dumbo
        dumbo.run(mapper, reducer)
