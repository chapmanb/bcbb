"""Calculate top cross species hits using BLAST.

Uses best e-value as a threshold to identify best cross-species hits in a number
of organism databases.
"""
import os
import subprocess
import contextlib
import tempfile
import xml.parsers.expat

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

def get_org_dbs(db_dir, target_org):
    """Retrieve references to fasta and BLAST databases for included organisms.
    """
    fasta_ref = None
    org_names = []
    db_refs = []
    with open(os.path.join(db_dir, "organism_dbs.txt")) as org_handle:
        for line in org_handle:
            org, db = line.rstrip("\r\n").split("\t")
            if db:
                if org == target_org:
                    assert fasta_ref is None
                    fasta_ref = os.path.join(db_dir, "%s.fa" % db)
                org_names.append(org)
                db_refs.append(os.path.join(db_dir, db))
    assert fasta_ref is not None, "Did not find base organism database"
    return fasta_ref, org_names, db_refs

def blast_top_hits(key, rec, db_refs, tmp_dir):
    """BLAST a fasta record against multiple DBs, returning top IDs and scores.
    """
    cur_id = _normalize_id(key)
    id_info = []
    score_info = []
    with _tmpfile(prefix="in", dir=tmp_dir) as ref_in:
        with open(ref_in, "w") as out_handle:
            out_handle.write(rec)
        for xref_db in db_refs:
            with _tmpfile(prefix="out", dir=tmp_dir) as blast_out:
                out_id, out_eval = _compare_by_blast(ref_in, xref_db, blast_out)
                id_info.append(out_id)
                score_info.append(out_eval)
    return cur_id, id_info, score_info

def _compare_by_blast(input_ref, xref_db, blast_out):
    """Compare all genes in an input file to the output database.
    """
    cl = NcbiblastpCommandline(query=input_ref, db=xref_db, out=blast_out,
            outfmt=5, num_descriptions=1, num_alignments=0)
    subprocess.check_call(str(cl).split())
    with open(blast_out) as blast_handle:
        try:
            rec = NCBIXML.read(blast_handle)
        except xml.parsers.expat.ExpatError:
            rec = None
        if rec and len(rec.descriptions) > 0:
            id_info = _normalize_id(rec.descriptions[0].title.split()[1])
            return id_info, rec.descriptions[0].bits
        else:
            return "", 0

def _normalize_id(id_info):
    if id_info.startswith("gi|"):
        parts = [p for p in id_info.split("|") if p]
        id_info = parts[-1]
    return id_info

@contextlib.contextmanager
def _tmpfile(*args, **kwargs):
    """Make a tempfile, safely cleaning up file descriptors on completion.
    """
    (fd, fname) = tempfile.mkstemp(*args, **kwargs)
    try:
        yield fname
    finally:
        os.close(fd)
        if os.path.exists(fname):
            os.remove(fname)
