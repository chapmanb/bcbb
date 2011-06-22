"""Converts Illumina SampleSheet CSV files to the run_info.yaml input file.

This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.
"""
import os
import sys
import csv
import itertools
import difflib
import glob
from Bio.Seq import Seq

import yaml

from bcbio.solexa.flowcell import (get_flowcell_info)
from bcbio import utils

def _organize_lanes(info_iter, barcode_ids):
    """Organize flat lane information into nested YAML structure.
    """
    all_lanes = []
    for (fcid, lane, sampleref), info in itertools.groupby(info_iter, lambda x: (x[0], x[1], x[1])):
        info = list(info)
        cur_lane = dict(flowcell_id=fcid, lane=lane, genome_build=info[0][3], analysis="Standard")
        
        if not _has_barcode(info):
            cur_lane["description"] = info[0][1]
        else: # barcoded sample
            cur_lane["description"] = "Lane %s, %s" % (lane, cur_lane["description"])
            multiplex = []
            for (_, _, sample_id, _, bc_seq) in info:
                bc_type, bc_id = barcode_ids[bc_seq]
                multiplex.append(dict(barcode_type=bc_type,
                                      barcode_id=bc_id,
                                      sequence=bc_seq,
                                      name=sample_id))
            cur_lane["multiplex"] = multiplex
        all_lanes.append(cur_lane)
    return all_lanes

def _has_barcode(sample):
    if sample[0][4]:
        return True
    else:
        raise "No barcode present on samplesheet sample %s !" % sample

def _generate_barcode_ids(info_iter):
    """Create unique barcode IDs assigned to sequences
    """
    bc_type = "SampleSheet"
    barcodes = list(set([x[-1] for x in info_iter]))
    barcodes.sort()
    barcode_ids = {}
    for i, bc in enumerate(barcodes):
        barcode_ids[bc] = (bc_type, i+1)
    return barcode_ids

def _read_input_csv(in_file):
    """Parse useful details from SampleSheet CSV file.
    """
    with open(in_file, "rU") as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for line in reader:
            if line: # empty lines
                (fc_id, lane, sample_id, genome, barcode) = line[:5]
                yield fc_id, lane, sample_id, genome, barcode

def _get_flowcell_id(in_file, require_single=True):
    """Retrieve the unique flowcell id represented in the SampleSheet.
    """
    fc_ids = set([x[0] for x in _read_input_csv(in_file)])
    if require_single and len(fc_ids) > 1:
        raise ValueError("There are several FCIDs in the same samplesheet file: %s" % in_file)
    else:
        return fc_ids

def _check_illumina_idx(sample_id, bc_seq):
    """ Sanity checks for barcodes: Makes sure "SampleID" matches the
        actual illumina sequence on the "Index" samplesheet column
    """

    official_indexes = {
        'index1': 'CGTGAT',
        'index2': 'ACATCG',
        'index3': 'GCCTAA',
        'index4': 'TGGTCA',
        'index5': 'CACTGT',
        'index6': 'ATTGGC',
        'index7': 'GATCTG',
        'index8': 'TCAAGT',
        'index9': 'CTGATC',
        'index10': 'AAGCTA',
        'index11': 'GTAGCC',
        'index12': 'TACAAG'
    }

    sample_idx = sample_id.split("_")
    if not len(sample_idx) == 1:
        sample_idx = sample_idx[-1]

    assert sample_idx in official_indexes.keys(), "SampleID %s does not conform *_indexN format" % sample

    # Official (non-custom) barcode ends in "A" ?
    if len(bc_seq) == 6:
        raise ValueError("Barcode %s needs a final A base" % bc_seq)
    if len(bc_seq) == 7:
        assert bc_seq[-1] == "A", "Barcode %s does not have a final A base, demultiplexing tainted ?" % bc_seq

    # We'll check for barcodes and its reverse complements
    official_idx = official_indexes[sample_idx]
    official_idx_rc = str(Seq(official_idx).reverse_complement())

    assert sample_idx in official_indexes.keys(), "Found SampleID %. Does not match any official illumina barcode."
    assert official_idx == bc_seq or official_idx_rc == bc_seq \
           or official_idx+"A" == bc_seq or official_idx_rc+"A" == bc_seq, \
           "\nOfficial illumina %s corresponds to %s or %s \
            \nSamplesheet reads %s corresponds to %s" % (sample_idx, official_idx, official_idx_rc,
                                                         sample_idx, bc_seq)

def csv2yaml(in_file, out_file=None):
    """Convert a CSV SampleSheet to YAML run_info format.
    """
    if out_file is None:
        out_file = "%s.yaml" % os.path.splitext(in_file)[0]

    samplesheet = _read_input_csv(in_file)
    sanity_checks(samplesheet)

    barcode_ids = _generate_barcode_ids(samplesheet)
    lanes = _organize_lanes(samplesheet, barcode_ids)
    with open(out_file, "w") as out_handle:
        out_handle.write(yaml.dump(lanes, default_flow_style=False))
    return out_file


def sanity_checks(samplesheet):
    for (_, _, sample_id, _, bc_seq) in samplesheet:
        _check_illumina_idx(sample_id, bc_seq)

def run_has_samplesheet(fc_dir, config, require_single=True):
    """Checks if there's a suitable SampleSheet.csv present for the run
    """
    fc_name, _ = get_flowcell_info(fc_dir)
    sheet_dirs = config.get("samplesheet_directories", [])
    fcid_sheet = {}
    for ss_dir in (s for s in sheet_dirs if os.path.exists(s)):
        with utils.chdir(ss_dir):
            for ss in glob.glob("*.csv"):
                fc_ids = _get_flowcell_id(ss, require_single)
                for fcid in fc_ids:
                    if fcid:
                        fcid_sheet[fcid] = os.path.join(ss_dir, ss)
    # difflib handles human errors while entering data on the SampleSheet.
    # Only one best candidate is returned (if any). 0.85 cutoff allows for
    # maximum of 2 mismatches in fcid

    potential_fcids = difflib.get_close_matches(fc_name, fcid_sheet.keys(), 1, 0.85)
    if len(potential_fcids) > 0 and fcid_sheet.has_key(potential_fcids[0]):
        return fcid_sheet[potential_fcids[0]]
    else:
        return None
