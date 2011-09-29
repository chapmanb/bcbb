"""Top level driver functionality for processing a sequencing lane.
"""
import os
import copy
import csv
import glob

from operator import itemgetter
from bcbio.utils import UnicodeWriter
from bcbio.pipeline import log
from bcbio.pipeline.fastq import get_fastq_files, get_multiplex_items
from bcbio.pipeline.demultiplex import split_by_barcode
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.solexa.flowcell import get_flowcell_info

def process_lane(info, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    config = _update_config_w_custom(config, info)

    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)

    log.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))

    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"], info, fc_name)
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    lane_items = []
    for mname, msample, fastq1, fastq2 in split_by_barcode(full_fastq1,
            full_fastq2, multiplex, lane_name, dirs, config):
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample,
                           dirs, config))
    
    _write_demultiplex_counts(lane_name,fc_name,fc_date,dirs["work"],info)
    
    return lane_items
    
# Append the demultiplexing results for this lane to the report file
def _write_demultiplex_counts(lane_name, fc_name, fc_date, workdir, info):
    
    multiplex = info.get("multiplex", None)
    if not multiplex or len(multiplex) == 0:
        return
    
    lane = info['lane']
    description = info.get('description','Lane '+str(info['lane']))
    
    metrics_file = os.path.join(workdir, "%s_barcode" % lane_name, "%s_bc.metrics" % lane_name)
    dmplx_report_file = os.path.join(workdir, "%s_%s_demultiplexed_read_counts.txt" % (fc_date,fc_name))
    
    if not os.path.exists(metrics_file):
        log.info("No barcode metrics file could be found for lane %s (expected '%s')" % (lane_name,metrics_file))
        return
    
    # Get the name and sequence corresponding to each barcode index in the multiplex section for this lane
    bc_meta_data = dict()
    for bc in multiplex:
        bc_id = bc.get("barcode_id",None)
        if bc_id:
            bc_meta_data[str(bc_id)] = [bc.get('name','N/A'),bc.get('sequence','N/A'),bc.get('barcode_type','N/A')]
            
    # Parse the demultiplexed barcode counts in the metrics file
    bc_metrics = []
    with open(metrics_file,"rb") as mfr:
        csvr = csv.reader(mfr,dialect='excel-tab')
        for row in csvr:
            bc_metrics.append([row[0],int(row[1])])
    
    # Sort the metrics according to number of reads
    bc_metrics = sorted(bc_metrics,key=itemgetter(1),reverse=True)
    
    # Join the metrics and the meta data into a list
    dmplx = []
    rank = 1
    for (bc_id,bc_count) in bc_metrics:
        
        # Get the corresponding meta data as specified in the multiplex section of run_info.yaml
        meta = bc_meta_data.get(str(bc_id),None)
        comment = ""
        
        # If meta data was not specified, check that the number of reads for this barcode isn't greater than the number for any of the expected barcodes
        if not meta:
            meta = ["N/A","N/A","N/A"]
            if rank <= len(multiplex):
                comment = "not specified in samplesheet, yet more reads than some of the expected barcodes"
        # Else, delete the entry from the barcode dictionary
        else:
            del bc_meta_data[str(bc_id)]
            if rank > len(multiplex):
                pass
    
        rank += 1
        dmplx.append([fc_date,fc_name,lane,description,bc_id,meta[0],meta[1],meta[2],bc_count,comment])
    
    # Warn if no demultiplex metrics were found for a barcode entry in run_info.yaml
    for bc_id, meta in bc_meta_data.items():
        dmplx.append([fc_date,fc_name,lane,description,bc_id,meta[0],meta[1],meta[2],0,"no reads after demultiplexing"])
    
    # Append the list of results to the report file 
    with open(dmplx_report_file,"ab") as mfw:
        csvw = UnicodeWriter(mfw, dialect='excel-tab')
        csvw.writerows(dmplx)    

def process_alignment(fastq1, fastq2, genome_build, lane_name, sample, dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    out_bam = ""
    if os.path.exists(fastq1) and aligner:
        log.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        out_bam = align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                                    lane_name, sample, dirs, config)
    return [{"sample": sample, "fastq": [fastq1, fastq2], "out_bam": out_bam,
            "lane": lane_name}]

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).iteritems():
        config["algorithm"][key] = val
    return config

def make_lane_items(info, fc_date, fc_name, dirs, config):
    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", "")
    log.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    lane_items = []
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
        mitems = get_multiplex_items(multiplex, info['lane'], dirs['fc_dir'], fc_name, fc_date)
        for fastq1, fastq2, mlane_name, msample in mitems:
            lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample, dirs, config))
    else:
        # TODO: Not multiplex: what to do?
        pass
    return lane_items

def get_flowcell_id(run_info, fc_dir, check_bc=True, glob_ext="_fastq.txt"):
    lane = None
    for info in run_info:
        lane = info.get("lane", "")
    if check_bc:
        glob_str = "%s_*_barcode/*%s" % (lane, glob_ext)
    else:
        glob_str = "%s_*%s" % (lane, glob_ext)
    files = glob.glob(os.path.join(fc_dir, glob_str))
    try:
        (name, date) = get_flowcell_info(os.path.basename(files[0]))
    except:
        raise StandardError("No flowcell information found in " + str(fc_dir))
    return name, date
