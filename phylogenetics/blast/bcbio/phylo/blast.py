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
    id_info = [cur_id]
    score_info = [cur_id]
    for xref_db in db_refs:
        with _blast_filenames(tmp_dir) as (in_file, out_file):
            with open(in_file, "w") as out_handle:
                out_handle.write(rec)
            out_id, out_eval = _compare_by_blast(in_file, xref_db, out_file)
            id_info.append(out_id)
            score_info.append(out_eval)
    return id_info, score_info

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
def _blast_filenames(tmp_dir):
    """Create files needed for blast and cleanup on closing.
    """
    (fd1, in_file) = tempfile.mkstemp(prefix="in", dir=tmp_dir)
    (fd2, out_file) = tempfile.mkstemp(prefix="out", dir=tmp_dir)
    try:
        yield (in_file, out_file)
    finally:
        os.close(fd1)
        os.close(fd2)
        for fname in [in_file, out_file]:
            if os.path.exists(fname):
                os.remove(fname)
