"""Retrieve GEO data for an experiment, classifying groups by expression data.
"""
import sys
import os
import csv
import collections
import json
import cPickle

from Bio import Entrez
import rpy2.robjects as robjects

def main():
    organism = "Mus musculus"
    cell_types = ["proB", "ProB", "pro-B"]
    email = "chapmanb@50mail.com"
    save_dir = os.getcwd()
    exp_data = get_geo_data(organism, cell_types, email, save_dir,
        _is_wild_type)

def _is_wild_type(result):
    """Check if a sample is wild type from the title.
    """
    return result.samples[0][0].startswith("WT")

def get_geo_data(organism, cell_types, email, save_dir, is_desired_result):
    save_file = os.path.join(save_dir, "%s-results.pkl" % cell_types[0])
    if not os.path.exists(save_file):
        results = cell_type_gsms(organism, cell_types, email)
        for result in results:
            if is_desired_result(result):
                with open(save_file, "w") as out_handle:
                    cPickle.dump(result, out_handle)
                break

    with open(save_file) as save_handle:
        result = cPickle.load(save_handle)
    print result
    exp = result.get_expression(save_dir)
    for gsm_id, exp_info in exp.items():
        print gsm_id, exp_info.items()[:5]
    return exp

class GEOResult:
    """Represent a GEO summary with associated samples, getting expression data.
    """
    def __init__(self, summary, samples):
        self.summary = summary
        self.samples = samples

    def __str__(self):
        out = "- %s\n" % self.summary
        for title, accession in self.samples:
            out += " %s %s\n" % (title, accession)
        return out

    def get_expression(self, save_dir):
        """Retrieve microarray results for our samples mapped to transcript IDs
        """
        results = dict()
        for (title, gsm_id) in self.samples:
            tx_to_exp = self.get_gsm_tx_values(gsm_id, save_dir)
            results[title] = tx_to_exp
        return results

    def get_gsm_tx_values(self, gsm_id, save_dir):
        """Retrieve a map of transcripts to expression from a GEO GSM file.
        """
        gsm_meta_file = os.path.join(save_dir, "%s-meta.txt" % gsm_id)
        gsm_table_file = os.path.join(save_dir, "%s-table.txt" % gsm_id)
        if (not os.path.exists(gsm_meta_file) or 
                not os.path.exists(gsm_table_file)):
            self._write_gsm_map(gsm_id, gsm_meta_file, gsm_table_file)

        with open(gsm_meta_file) as in_handle:
            gsm_meta = json.load(in_handle)
        id_to_tx = self.get_gpl_map(gsm_meta['platform_id'], save_dir)
        tx_to_vals = collections.defaultdict(list)
        with open(gsm_table_file) as in_handle:
            reader = csv.reader(in_handle, dialect='excel-tab')
            reader.next() # header
            for probe_id, probe_val in reader:
                for tx_id in id_to_tx.get(probe_id, []):
                    tx_to_vals[tx_id].append(float(probe_val))
        return tx_to_vals

    def _write_gsm_map(self, gsm_id, meta_file, table_file):
        """Retrieve GEO expression values using Bioconductor, saving to a table.
        """
        robjects.r.assign("gsm.id", gsm_id)
        robjects.r.assign("table.file", table_file)
        robjects.r.assign("meta.file", meta_file)
        robjects.r('''
          library(GEOquery)
          library(rjson)
          gsm <- getGEO(gsm.id)
          write.table(Table(gsm), file = table.file, sep = "\t", row.names = FALSE,
                      col.names = TRUE)
          cat(toJSON(Meta(gsm)), file = meta.file)
        ''')

    def get_gpl_map(self, gpl_id, save_dir):
        """Retrieve a map of IDs to transcript information from a GEO GPL file.
        """
        gpl_file = os.path.join(save_dir, "%s-map.txt" % gpl_id)
        if not os.path.exists(gpl_file):
            self._write_gpl_map(gpl_id, gpl_file)
        gpl_map = collections.defaultdict(list)
        with open(gpl_file) as in_handle:
            reader = csv.reader(in_handle, dialect='excel-tab')
            reader.next() # header
            for probe_id, tx_id_str in reader:
                for tx_id in tx_id_str.split(' /// '):
                    if tx_id:
                        gpl_map[probe_id].append(tx_id)
        return dict(gpl_map)

    def _write_gpl_map(self, gpl_id, gpl_file):
        """Retrieve GEO platform data using R and save to a table.
        """
        robjects.r.assign("gpl.id", gpl_id)
        robjects.r.assign("gpl.file", gpl_file)
        robjects.r('''
          library(GEOquery)
          gpl <- getGEO(gpl.id)
          gpl.map <- subset(Table(gpl), select=c("ID", "RefSeq.Transcript.ID"))
          write.table(gpl.map, file = gpl.file, sep = "\t", row.names = FALSE,
                      col.names = TRUE)
        ''')

def cell_type_gsms(organism, cell_types, email):
    """Use Entrez to retrieve GEO entries for an organism and cell type.
    """
    Entrez.email = email
    search_term = "%s[ORGN] %s" % (organism, " OR ".join(cell_types))
    print "Searching GEO and retrieving results: %s" % search_term
    
    hits = []
    handle = Entrez.esearch(db="gds", term=search_term)
    results = Entrez.read(handle)
    for geo_id in results['IdList']:
        handle = Entrez.esummary(db="gds", id=geo_id)
        summary = Entrez.read(handle)
        samples = []
        for sample in summary[0]['Samples']:
            for cell_type in cell_types:
                if sample['Title'].find(cell_type) >= 0:
                    samples.append((sample['Title'], sample['Accession']))
                    break
        if len(samples) > 0:
            hits.append(GEOResult(summary[0]['summary'], samples))
    return hits

if __name__ == "__main__":
    main(*sys.argv[1:])
