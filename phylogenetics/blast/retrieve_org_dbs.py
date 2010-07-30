#!/usr/bin/env python
"""Retrieve full genome databases, preparing them for BLAST analysis.

Usage:
    retrieve_org_dbs.py <YAML config file>

Requires:
    - NCBI's blast+ -- for preparing the organism databases
      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
"""
import os
import sys
import csv
import ftplib
import subprocess
import contextlib

import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    retriever = EnsemblFtpRetrieval()
    organisms = read_org_list(config['org_file'])
    db_dir = config['db_dir']
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    org_files = []
    for org in organisms:
        db_file = retriever.retrieve_db(org, db_dir)
        org_files.append((org, db_file))
    with open(os.path.join(db_dir, "organism_dbs.txt"), "w") as out_handle:
        for org, fname in org_files:
            out_handle.write("%s\t%s\n" % (org, fname))

def read_org_list(in_file):
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        orgs = [r[-1] for r in reader]
    return orgs

class EnsemblFtpRetrieval:
    """Handle obtaining a reference genome from Ensembl
    """
    def __init__(self):
        self._main_ftp = "ftp://ftp.ensembl.org/pub/current_fasta/"
        self._genome_ftp = "ftp://ftp.ensemblgenomes.org/pub/%s/current/fasta/"
        self._genome_dbs = ["bacteria", "protists", "metazoa", "fungi",
                "plants"]
        urls = [self._genome_ftp % d for d in self._genome_dbs] + \
               [self._main_ftp]
        self._org_to_urls = dict()
        for url in urls:
            orgs = self._files_at_url(url)
            for org in orgs:
                self._org_to_urls[org] = url

    def _files_at_url(self, url):
        """Add organisms available at the provided FTP url.
        """
        parts = url.replace("ftp://", "").split("/")
        ftp = ftplib.FTP(parts[0])
        ftp.login()
        orgs = ftp.nlst("/".join(parts[1:]))
        return [o.split("/")[-1] for o in orgs]

    def retrieve_db(self, organism, db_dir):
        ftp_url = self._get_ftp_url(organism)
        if ftp_url is None:
            return ""
        file_name = ftp_url.split("/")[-1]
        final_file = os.path.join(db_dir, file_name.replace(".gz", ""))
        db_name = os.path.splitext(os.path.basename(final_file))[0]
        if not os.path.exists(final_file):
            with _chdir(db_dir):
                cl = ["wget", ftp_url]
                subprocess.check_call(cl)
                cl = ["gunzip", file_name]
                subprocess.check_call(cl)
                cl = ["makeblastdb", "-in", os.path.basename(final_file),
                      "-dbtype", "prot",
                      "-out", db_name,
                      "-title", organism]
                subprocess.check_call(cl)
        return db_name

    def _get_ftp_url(self, organism):
        """Retrieve the protein database link for a given organism.
        """
        ftp_url = None
        org_parts = organism.split()
        for check_org in [organism.replace(" ", "_").lower(),
                "_".join([org_parts[0][0], org_parts[1]]).lower()]:
            try:
                ftp_url = self._org_to_urls[check_org]
                break
            except KeyError:
                pass
        if ftp_url:
            ftp_url = ftp_url + check_org + "/pep/"
            files = self._files_at_url(ftp_url)
            for f in files:
                if f.endswith("pep.all.fa.gz"):
                    ftp_url = ftp_url + f
                    break
        return ftp_url

@contextlib.contextmanager
def _chdir(new_dir):
    orig_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(orig_dir)


if __name__ == "__main__":
    main(*sys.argv[1:])
