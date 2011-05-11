#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file> <flow cell dir>
                                  [<YAML run information>]

The optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml'

Workflow:
    - Retrieve details on a run.
    - Align fastq files to reference genome.
    - Perform secondary analyses like SNP calling.
    - Generate summary report.
"""
import os
import sys
import json
import contextlib
import subprocess
import glob
import copy
import csv
import copy
import shutil
import collections
from optparse import OptionParser
import xml.etree.ElementTree as ET
import StringIO

import yaml
import logbook
from Bio import SeqIO

from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir)
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.broad.metrics import PicardMetricsParser
from bcbio import utils
from bcbio.broad import BroadRunner
from bcbio.log import create_log_handler

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(config_file, fc_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, LOG_NAME)
    
    with log_handler.applicationbound():
        run_main(config, config_file, fc_dir, run_info_yaml)

def run_main(config, config_file, fc_dir, run_info_yaml):
    work_dir = os.getcwd()
    fc_name, fc_date = get_flowcell_info(fc_dir)

    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        run_info = dict(details=run_details, run_id="")
    else:
        log.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_info = galaxy_api.run_details(fc_name, fc_date)
    fastq_dir = get_fastq_dir(fc_dir)
    run_items = _add_multiplex_across_lanes(run_info["details"], fastq_dir, fc_name)
    align_dir = os.path.join(work_dir, "alignments")

    # process each flowcell lane
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_lane,
                       ((i, fastq_dir, fc_name, fc_date, align_dir, config, config_file)
                        for i in run_items)):
            pass
    # process samples, potentially multiplexed across multiple lanes
    sample_files, sample_fastq, sample_info = organize_samples(align_dir,
            fastq_dir, work_dir, fc_name, fc_date, run_items)
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_sample, ((name, sample_fastq[name], sample_info[name],
                                         bam_files, work_dir, config, config_file)
                                        for name, bam_files in sample_files)):
            pass
    write_metrics(run_info, work_dir, fc_dir, fc_name, fc_date, fastq_dir)

@utils.map_wrap
def process_lane(info, fastq_dir, fc_name, fc_date, align_dir, config,
        config_file):
    """Do alignments for a lane, potentially splitting based on barcodes.
    """
    config = _update_config_w_custom(config, info)
    
    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)
    
    log.info("Processing sample %s on lane %s with reference genome %s by researcher %s. Using %s analysis preset" \
             % (sample_name, info["lane"], genome_build, \
                info.get("researcher", "unknown"), info.get("analysis", "")))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
    
    fastq_dir, galaxy_dir = _get_full_paths(fastq_dir, config, config_file)
    align_ref, sam_ref = get_genome_ref(genome_build,
            config["algorithm"]["aligner"], galaxy_dir)
    full_fastq1, full_fastq2 = get_fastq_files(fastq_dir, info['lane'], fc_name)
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    for mname, msample, fastq1, fastq2 in split_by_barcode(full_fastq1,
            full_fastq2, multiplex, lane_name, config):
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        if os.path.exists(fastq1) and config["algorithm"]["aligner"]:
            do_alignment(fastq1, fastq2, align_ref, sam_ref, mlane_name,
                    msample, align_dir, config, config_file)

@utils.map_wrap
def process_sample(sample_name, fastq_files, info, bam_files, work_dir,
        config, config_file):
    """Finalize processing for a sample, potentially multiplexed.
    """
    config = _update_config_w_custom(config, info)
    
    genome_build = info.get("genome_build", None)
    (_, galaxy_dir) = _get_full_paths("", config, config_file)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  galaxy_dir)
    fastq1, fastq2 = _combine_fastq_files(fastq_files, work_dir)
    log.info("Combining and preparing wig file %s" % str(sample_name))
    sort_bam = merge_bam_files(bam_files, work_dir, config)
    bam_to_wig(sort_bam, config, config_file)
    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        dbsnp_file = get_dbsnp_file(config, sam_ref)
        gatk_bam = recalibrate_quality(sort_bam, sam_ref,
                dbsnp_file, config_file)
        log.info("Analyzing recalibration %s" % str(sample_name))
        analyze_recalibration(gatk_bam, fastq1, fastq2)
        if config["algorithm"]["snpcall"]:
            log.info("SNP genotyping %s with GATK" % str(sample_name))
            vrn_file = run_genotyper(gatk_bam, sam_ref, dbsnp_file, config_file)
            eval_genotyper(vrn_file, sam_ref, dbsnp_file, config)
            log.info("Calculating variation effects for %s" % str(sample_name))
            variation_effects(vrn_file, genome_build, sam_ref, config)
    if sam_ref is not None:
        print sample_name, "Generating summary files"
        generate_align_summary(sort_bam, fastq2 is not None, sam_ref,
                config, sample_name, config_file)

def _combine_fastq_files(in_files, work_dir):
    if len(in_files) == 1:
        return in_files[0]
    else:
        cur1, cur2 = in_files[0]
        out1 = os.path.join(work_dir, os.path.basename(cur1))
        out2 = os.path.join(work_dir, os.path.basename(cur2)) if cur2 else None
        if not os.path.exists(out1):
            with open(out1, "a") as out_handle:
                for (cur1, _) in in_files:
                    with open(cur1) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        if out2 and not os.path.exists(out2):
            with open(out2, "a") as out_handle:
                for (_, cur2) in in_files:
                    with open(cur2) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        return out1, out2

def organize_samples(align_dir, fastq_dir, work_dir, fc_name, fc_date, run_items):
    """Organize BAM output files by sample name, handling multiplexing.
    """
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    for lane_info in run_items:
        multiplex = lane_info.get("multiplex", None)
        if multiplex:
            mfastq_dir = os.path.join(work_dir, "%s_%s_%s_barcode" %
                    (lane_info["lane"], fc_date, fc_name))
            for multi in multiplex:
                name = (lane_info.get("name", ""), lane_info["description"],
                        multi["name"])
                fname = os.path.join(align_dir, "%s_%s_%s_%s-sort.bam" %
                    (lane_info["lane"], fc_date, fc_name, multi["barcode_id"]))
                if os.path.exists(fname):
                    bams_by_sample[name].append(fname)
                    sample_info[name] = lane_info
                    fastq_by_sample[name].append(get_fastq_files(mfastq_dir,
                        lane_info["lane"], fc_name, multi["barcode_id"]))
        else:
            name = (lane_info.get("name", ""), lane_info["description"])
            fname = os.path.join(align_dir, "%s_%s_%s-sort.bam" %
                    (lane_info["lane"], fc_date, fc_name))
            if os.path.exists(fname):
                bams_by_sample[name].append(fname)
                sample_info[name] = lane_info
                fastq_by_sample[name].append(get_fastq_files(fastq_dir,
                    lane_info["lane"], fc_name))
    return sorted(bams_by_sample.items()), dict(fastq_by_sample), sample_info

def _get_fastq_size(item, fastq_dir, fc_name):
    """Retrieve the size of reads from the first flowcell sequence.
    """
    (fastq1, _) = get_fastq_files(fastq_dir, item['lane'], fc_name)
    with open(fastq1) as in_handle:
        try:
            rec = SeqIO.parse(in_handle, "fastq").next()
            size = len(rec.seq)
        except StopIteration:
            log.warn("Found a zero-sized fastq file for lane %s" % item['lane'])
            size = 0
    return size

def _add_multiplex_across_lanes(run_items, fastq_dir, fc_name):
    """Add multiplex information to control and non-multiplexed lanes.

    Illumina runs include barcode reads for non-multiplex lanes, and the
    control, when run on a multiplexed flow cell. This checks for this
    situation and adds details to trim off the extra bases.
    """
    fastq_dir = _add_full_path(fastq_dir)
    # determine if we have multiplexes and collect expected size
    fastq_sizes = []
    tag_sizes = []
    has_barcodes = False
    for item in run_items:
        if item.get("multiplex", None):
            has_barcodes = True
            tag_sizes.extend([len(b["sequence"]) for b in item["multiplex"]])
            fastq_sizes.append(_get_fastq_size(item, fastq_dir, fc_name))
    if not has_barcodes: # nothing to worry about
        return run_items
    fastq_sizes = list(set(fastq_sizes))

    # discard 0 sizes to handle the case where lane(s) are empty or failed
    try:
        fastq_sizes.remove(0)
    except ValueError: pass
    
    tag_sizes = list(set(tag_sizes))
    final_items = []
    for item in run_items:
        if item.get("multiplex", None) is None:
            assert len(fastq_sizes) == 1, \
                   "Multi and non-multiplex reads with multiple sizes"
            expected_size = fastq_sizes[0]
            assert len(tag_sizes) == 1, \
                   "Expect identical tag size for a flowcell"
            tag_size = tag_sizes[0]
            this_size = _get_fastq_size(item, fastq_dir, fc_name)
            if this_size == expected_size:
                item["multiplex"] = [{"name" : item.get("name", item["description"]),
                                      "barcode_id": "trim",
                                      "sequence" : "N" * tag_size}]
            else:
                assert this_size == expected_size - tag_size, \
                       "Unexpected non-multiplex sequence"
        final_items.append(item)
    return final_items

def split_by_barcode(fastq1, fastq2, multiplex, base_name, config):
    """Split a fastq file into multiplex pieces using barcode details.
    """
    if not multiplex:
        return [("", "", fastq1, fastq2)]
    bc_dir = "%s_barcode" % base_name
    nomatch_file = "%s_unmatched_1_fastq.txt" % base_name
    metrics_file = "%s_bc.metrics" % base_name
    with utils.chdir(bc_dir):
        if not os.path.exists(nomatch_file) and not os.path.exists(metrics_file):
            tag_file = _make_tag_file(multiplex)
            cl = [config["program"]["barcode"], tag_file,
                  "%s_--b--_--r--_fastq.txt" % base_name,
                  fastq1]
            if fastq2:
                cl.append(fastq2)
            cl.append("--mismatch=%s" % config["algorithm"]["bc_mismatch"])
            cl.append("--metrics=%s" % metrics_file)
            if int(config["algorithm"]["bc_read"]) == 2:
                cl.append("--second")
            if int(config["algorithm"]["bc_position"]) == 5:
                cl.append("--five")
            subprocess.check_call(cl)
    out_files = []
    for info in multiplex:
        fq_fname = lambda x: os.path.join(bc_dir, "%s_%s_%s_fastq.txt" %
                             (base_name, info["barcode_id"], x))
        bc_file1 = fq_fname("1")
        bc_file2 = fq_fname("2") if fastq2 else None
        out_files.append((info["barcode_id"], info["name"], bc_file1, bc_file2))
    return out_files

def _make_tag_file(barcodes):
    tag_file = "%s-barcodes.cfg" % barcodes[0].get("barcode_type", "barcode")
    barcodes = _adjust_illumina_tags(barcodes)
    with open(tag_file, "w") as out_handle:
        for bc in barcodes:
            out_handle.write("%s %s\n" % (bc["barcode_id"], bc["sequence"]))
    return tag_file

def _adjust_illumina_tags(barcodes):
    """Handle additional trailing A in Illumina barocdes.

    Illumina barcodes are listed as 6bp sequences but have an additional
    A base when coming off on the sequencer. This checks for this case and
    adjusts the sequences appropriately if needed.
    """
    illumina_size = 7
    all_illumina = True
    need_a = False
    for bc in barcodes:
        if bc.get("barcode_type", "illumina").lower().find("illumina") == -1:
            all_illumina = False
        if (not bc["sequence"].upper().endswith("A") or
            len(bc["sequence"]) < illumina_size):
            need_a = True
    if all_illumina and need_a:
        new = []
        for bc in barcodes:
            new_bc = copy.deepcopy(bc)
            new_bc["sequence"] = "%sA" % new_bc["sequence"]
            new.append(new_bc)
        barcodes = new
    return barcodes

def do_alignment(fastq1, fastq2, align_ref, sam_ref, lane_name,
        sample_name, align_dir, config, config_file):
    """Align to the provided reference genome, returning an aligned SAM file.
    """
    aligner_to_use = config["algorithm"]["aligner"]
    utils.safe_makedir(align_dir)

    log.info("Aligning lane %s with %s aligner" % (lane_name, aligner_to_use))
    if aligner_to_use == "bowtie":
        sam_file = bowtie_to_sam(fastq1, fastq2, align_ref, lane_name,
                                 align_dir, config)
    elif aligner_to_use == "tophat":
        sam_file = tophat_align_to_sam(fastq1, fastq2, align_ref, lane_name,
                                       align_dir, config)
    elif aligner_to_use == "bwa":
        sam_file = bwa_align_to_sam(fastq1, fastq2, align_ref, lane_name,
                                    align_dir, config)
    elif aligner_to_use == "maq":
        sam_file = maq_align_to_sam(fastq1, fastq2, align_ref, lane_name,
                                    sample_name, align_dir, config_file)
    else:
        raise ValueError("Do not recognize aligner: %s" % aligner_to_use)
    log.info("Converting lane %s to sorted BAM file" % lane_name)
    sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, sample_name,
                    lane_name, config, config_file)

def bowtie_to_sam(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    """Before a standard or paired end alignment with bowtie.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file):
        cl = [config["program"]["bowtie"]]
        cl += ["-q", "--solexa1.3-quals",
               "-v", config["algorithm"]["max_errors"],
               "-k", 1,
               "-X", 1000, # matches bwa sampe default size
               "-M", 1,
               "--best",
               "--strata",
               "--sam",
               ref_file]
        if pair_file:
            cl += ["-1", fastq_file, "-2", pair_file]
        else:
            cl += [fastq_file]
        cl += [out_file]
        cl = [str(i) for i in cl]
        log.info("Running bowtie with cmdline: %s" % " ".join(cl))
        subprocess.check_call(cl)
    return out_file

def tophat_align_to_sam(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    out_dir = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_dir):
        cl = [config["program"]["tophat"]]
        cl += ["--solexa1.3-quals",
              "-p 8",
              "-r 45",
              ref_file]
        if pair_file:
            cl += [fastq_file, pair_file]
        else:
            cl += [fastq_file]
        cl += ["-o ", out_dir]
    log.info("Running tophat with cmdline: %s" % cl)
    child = subprocess.check_call(cl)
    return out_dir

def bwa_align_to_sam(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    """Perform a BWA alignment, generating a SAM file.
    """
    sai1_file = os.path.join(align_dir, "%s_1.sai" % out_base)
    sai2_file = (os.path.join(align_dir, "%s_2.sai" % out_base)
                 if pair_file else None)
    sam_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(sam_file):
        if not os.path.exists(sai1_file):
            _run_bwa_align(fastq_file, ref_file, sai1_file, config)
        if sai2_file and not os.path.exists(sai2_file):
            _run_bwa_align(pair_file, ref_file, sai2_file, config)
        align_type = "sampe" if sai2_file else "samse"
        sam_cl = [config["program"]["bwa"], align_type, ref_file, sai1_file]
        if sai2_file:
            sam_cl.append(sai2_file)
        sam_cl.append(fastq_file)
        if sai2_file:
            sam_cl.append(pair_file)
        with open(sam_file, "w") as out_handle:
            subprocess.check_call(sam_cl, stdout=out_handle)
    return sam_file

def maq_align_to_sam(fastq_file, pair_file, ref_file, out_base,
        sample_name, align_dir, config_file):
    """Produce a BAM output using Picard to do a maq alignment.
    """
    raise NotImplementedError("Need to update for alignment directory")
    bam_file = os.path.join(align_dir, "%s-maq-cal.bam" % out_base)
    cl = ["picard_maq_recalibrate.py", "--name=%s" % sample_name,
            config_file, out_base, ref_file, fastq_file]
    if pair_file:
        cl.append(pair_file)
    child = subprocess.Popen(cl)
    child.wait()
    return bam_file

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config["program"]["bwa"], "aln",
              "-n %s" % config["algorithm"]["max_errors"],
              "-k %s" % config["algorithm"]["max_errors"],
              ref_file, fastq_file]
    with open(out_file, "w") as out_handle:
        subprocess.check_call(aln_cl, stdout=out_handle)

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, sample_name,
                    lane_name, config, config_file):
    """Convert SAM file to merged and sorted BAM file.
    """
    lane = lane_name.split("_")[0]
    cl = ["picard_sam_to_bam.py", "--name=%s" % sample_name,
            "--rg=%s" % lane, "--pu=%s" % lane_name,
            config_file, sam_file, ref_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)
    utils.save_diskspace(sam_file, "SAM converted to BAM", config)

def merge_bam_files(bam_files, work_dir, config):
    """Merge multiple BAM files from a sample into a single BAM for processing.
    """
    out_file = os.path.join(work_dir, os.path.basename(bam_files[0]))
    if not os.path.exists(out_file):
        picard = BroadRunner(config["program"]["picard"],
                             max_memory=config["algorithm"].get("java_memory", ""))
        with utils.curdir_tmpdir() as tmp_dir:
            opts = [("OUTPUT", out_file),
                    ("SORT_ORDER", "coordinate"),
                    ("TMP_DIR", tmp_dir)]
            for b in bam_files:
                opts.append(("INPUT", b))
            picard.run("MergeSamFiles", opts)
    for b in bam_files:
        utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
    return out_file

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        subprocess.check_call(cl)
    return wig_file

def generate_align_summary(bam_file, is_paired, sam_ref, config,
        sample_name, config_file, do_sort=False):
    """Run alignment summarizing script to produce a pdf with align details.
    """
    sample_name = " : ".join(sample_name)
    cl = ["align_summary_report.py", "--name=%s" % sample_name,
            config["program"]["picard"], bam_file, sam_ref]
    if is_paired:
        cl.append("--paired")
    bait = config["algorithm"].get("hybrid_bait", "")
    target = config["algorithm"].get("hybrid_target", "")
    if bait and target:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        cl.append("--bait=%s" % os.path.join(base_dir, bait))
        cl.append("--target=%s" % os.path.join(base_dir, target))
    if do_sort:
        cl.append("--sort")
    cl.append("--config=%s" % config_file)
    subprocess.check_call(cl)

def recalibrate_quality(sort_bam_file, sam_ref, dbsnp_file, config_file):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    bam_file = sort_bam_file.replace("-sort.bam", ".bam")
    cl = ["picard_gatk_recalibrate.py", config_file, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_files = glob.glob("%s*-gatkrecal.bam" % os.path.splitext(sort_bam_file)[0])
    assert len(out_files) == 1, out_files
    return out_files[0]

def run_genotyper(bam_file, ref_file, dbsnp_file, config_file):
    """Perform SNP genotyping and analysis using GATK.
    """
    cl = ["gatk_genotyper.py", config_file, ref_file, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    base, ext = os.path.splitext(bam_file)
    return glob.glob("%s*snp-filter.vcf" % base)[0]

def eval_genotyper(vrn_file, ref_file, dbsnp_file, config):
    """Evaluate variant genotyping, producing a JSON metrics file with values.
    """
    metrics_file = "%s.eval_metrics" % vrn_file
    cl = ["gatk_variant_eval.py", config["program"].get("gatk", config["program"]["picard"]),
          vrn_file, ref_file, dbsnp_file]
    target = config["algorithm"].get("hybrid_target", "")
    if target:
        base_dir = os.path.dirname(os.path.dirname(ref_file))
        cl.append(os.path.join(base_dir, target))
    with open(metrics_file, "w") as out_handle:
        subprocess.check_call(cl, stdout=out_handle)

def variation_effects(vrn_file, genome_build, ref_file, config):
    """Calculate effects of variations, associating them with transcripts.
    """
    snp_eff_dir = config["program"]["snpEff"]
    snp_eff_jar = os.path.join(snp_eff_dir, "snpEff.jar")
    cl = ["variant_effects.py", snp_eff_jar, vrn_file, genome_build]
    target = config["algorithm"].get("hybrid_target", "")
    if target:
        base_dir = os.path.dirname(os.path.dirname(ref_file))
        cl.append(os.path.join(base_dir, target))
    subprocess.check_call(cl)

def get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file

def analyze_recalibration(recal_file, fastq1, fastq2):
    """Provide a pdf report of GATK recalibration of scores.
    """
    cl = ["analyze_quality_recal.py", recal_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)

def get_fastq_files(directory, lane, fc_name, bc_name=None):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if bc_name:
        glob_str = "%s_*%s_%s_*_fastq.txt" % (lane, fc_name, bc_name)
    else:
        glob_str = "%s_*%s*_fastq.txt" % (lane, fc_name)
    files = glob.glob(os.path.join(directory, glob_str))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

def _remap_to_maq(ref_file):
    """ToDo: Why is this needed for ?
    """
    base_dir = os.path.dirname(os.path.dirname(ref_file))
    name = os.path.basename(ref_file)
    for ext in ["fa", "fasta"]:
        test_file = os.path.join(base_dir, "maq", "%s.%s" % (name, ext))
        if os.path.exists(test_file):
            return test_file
    raise ValueError("Did not find maq file %s" % ref_file)

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    if not aligner or not genome_build:
        return (None, None)
    ref_files = dict(
            bowtie = "bowtie_indices.loc",
            bwa = "bwa_index.loc",
            samtools = "sam_fa_indices.loc",
            maq = "bowtie_indices.loc",
	        tophat = "bowtie_indices.loc")
    remap_fns = dict(
            maq = _remap_to_maq
            )
    out_info = []
    ref_dir = os.path.join(galaxy_base, "tool-data")
    for ref_get in [aligner, "samtools"]:
        ref_file = os.path.join(ref_dir, ref_files[ref_get])
        cur_ref = None
        with open(ref_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        cur_ref = parts[-1]
                        break
        if cur_ref is None:
            raise IndexError("Genome %s not found in %s" % (genome_build,
                ref_file))
        try:
            cur_ref = remap_fns[ref_get](cur_ref)
        except KeyError:
            pass
        out_info.append(_add_full_path(cur_ref, ref_dir))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)

# Output high level summary information for a sequencing run in YAML format
# that can be picked up and loaded into Galaxy.

def write_metrics(run_info, analysis_dir, fc_dir, fc_name, fc_date,
        fastq_dir):
    """Write an output YAML file containing high level sequencing metrics.
    """
    lane_stats, sample_stats, tab_metrics = summary_metrics(run_info,
            analysis_dir, fc_name, fc_date, fastq_dir)
    out_file = os.path.join(analysis_dir, "run_summary.yaml")
    with open(out_file, "w") as out_handle:
        metrics = dict(lanes=lane_stats, samples=sample_stats)
        yaml.dump(metrics, out_handle, default_flow_style=False)
    tab_out_file = os.path.join(fc_dir, "run_summary.tsv")
    with open(tab_out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect="excel-tab")
        for info in tab_metrics:
            writer.writerow(info)
    return out_file

def summary_metrics(run_info, analysis_dir, fc_name, fc_date, fastq_dir):
    """Reformat run and analysis statistics into a YAML-ready format.
    """
    tab_out = []
    lane_info = []
    sample_info = []
    for run in run_info["details"]:
        tab_out.append([run["lane"], run.get("researcher", ""),
            run.get("name", ""), run.get("description")])
        base_info = dict(
                researcher = run.get("researcher_id", ""),
                sample = run.get("sample_id", ""),
                lane = run["lane"],
                request = run_info["run_id"])
        cur_lane_info = copy.deepcopy(base_info)
        cur_lane_info["metrics"] = _bustard_stats(run["lane"], fastq_dir,
                fc_date)
        lane_info.append(cur_lane_info)
        for barcode in run.get("multiplex", [None]):
            cur_name = "%s_%s_%s" % (run["lane"], fc_date, fc_name)
            if barcode:
                cur_name = "%s_%s-" % (cur_name, barcode["barcode_id"])
            stats = _metrics_from_stats(_lane_stats(cur_name, analysis_dir))
            if stats:
                cur_run_info = copy.deepcopy(base_info)
                cur_run_info["metrics"] = stats
                cur_run_info["barcode_id"] = str(barcode["barcode_id"]) if barcode else ""
                cur_run_info["barcode_type"] = (str(barcode.get("barcode_type", ""))
                                                if barcode else "")
                sample_info.append(cur_run_info)
    return lane_info, sample_info, tab_out

def _metrics_from_stats(stats):
    """Remap Broad metrics names to our local names.
    """
    if stats:
        s_to_m = dict(
                AL_MEAN_READ_LENGTH = 'Read length',
                AL_TOTAL_READS = 'Reads',
                AL_PF_READS_ALIGNED = 'Aligned',
                DUP_READ_PAIR_DUPLICATES = 'Pair duplicates'
                )
        metrics = dict()
        for stat_name, metric_name in s_to_m.iteritems():
            metrics[metric_name] = stats[stat_name]
        return metrics

def _bustard_stats(lane_num, fastq_dir, fc_date):
    """Extract statistics about the flow cell from Bustard outputs.
    """
    sum_file = os.path.join(fastq_dir, os.pardir, "BustardSummary.xml")
    #sum_file = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls",
    #        "BustardSummary.xml")
    stats = dict()
    if os.path.exists(sum_file):
        with open(sum_file) as in_handle:
            results = ET.parse(in_handle).getroot().find("TileResultsByLane")
            for lane in results:
                if lane.find("laneNumber").text == str(lane_num):
                    stats = _collect_cluster_stats(lane)
    read_stats = _calc_fastq_stats(fastq_dir, lane_num, fc_date)
    stats.update(read_stats)
    return stats

def _calc_fastq_stats(fastq_dir, lane_num, fc_date):
    """Grab read length from fastq; could provide distribution if non-equal.
    """
    stats = dict()
    fastq_files = glob.glob(os.path.join(fastq_dir, "%s_%s*" % (lane_num,
        fc_date)))
    if len(fastq_files) > 0:
        fastq_file = sorted(fastq_files)[-1]
        with open(fastq_file) as in_handle:
            line = in_handle.readline()
            stats["Read length"] = len(in_handle.readline().strip())
    return stats

def _collect_cluster_stats(lane):
    """Retrieve total counts on cluster statistics.
    """
    stats = {"Clusters" : 0, "Clusters passed": 0}
    for tile in lane.find("Read").findall("Tile"):
        stats["Clusters"] += int(tile.find("clusterCountRaw").text)
        stats["Clusters passed"] += int(tile.find("clusterCountPF").text)
    return stats

def _lane_stats(cur_name, work_dir):
    """Parse metrics information from files in the working directory.
    """
    parser = PicardMetricsParser()
    metrics_files = glob.glob(os.path.join(work_dir, "%s*metrics" % cur_name))
    metrics = parser.extract_metrics(metrics_files)
    return metrics

def _add_full_path(dirname, basedir=None):
    if basedir is None:
        basedir = os.getcwd()
    if not dirname.startswith("/"):
        dirname = os.path.join(basedir, dirname)
    return dirname

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = _add_full_path(fastq_dir)
    config_dir = _add_full_path(os.path.dirname(config_file))
    galaxy_config_file = _add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file)

# Utility functions

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    return config

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)
