"""Provide SNP, indel calling and variation analysis using GATK genotyping tools.

Genotyping:

http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3
http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper

Variant Evaluation:

http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
"""
import os
import copy
import itertools
import collections

from bcbio import broad
from bcbio.log import logger
from bcbio.utils import file_exists, replace_suffix, append_stem
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed.split import parallel_split_combine
from bcbio.pipeline.shared import (process_bam_by_chromosome, configured_ref_file,
                                   subset_variant_regions)
from bcbio.variation.realign import has_aligned_reads
import sh

# ## GATK Genotype calling

def _shared_gatk_call_prep(align_bam, ref_file, config, dbsnp, region, out_file):
    """Shared preparation work for GATK variant calling.
    """
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_fn("picard_index_ref", ref_file)
    broad_runner.run_fn("picard_index", align_bam)
    coverage_depth = config["algorithm"].get("coverage_depth", "high").lower()
    variant_regions = config["algorithm"].get("variant_regions", None)
    confidence = "4.0" if coverage_depth in ["low"] else "30.0"
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bam)[0]
    region = subset_variant_regions(variant_regions, region, out_file)

    params = ["-I", align_bam,
              "-R", ref_file,
              "--annotation", "QualByDepth",
              "--annotation", "HaplotypeScore",
              "--annotation", "MappingQualityRankSumTest",
              "--annotation", "ReadPosRankSumTest",
              "--annotation", "FisherStrand",
              "--annotation", "RMSMappingQuality",
              "--annotation", "DepthOfCoverage",
              "--standard_min_confidence_threshold_for_calling", confidence,
              "--standard_min_confidence_threshold_for_emitting", confidence,
              ]
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    if region:
        params += ["-L", region, "--interval_set_rule", "INTERSECTION"]
    return broad_runner, params, out_file

def unified_genotyper(align_bam, ref_file, config, dbsnp=None,
                       region=None, out_file=None):
    """Perform SNP genotyping on the given alignment file.
    """
    broad_runner, params, out_file = \
        _shared_gatk_call_prep(align_bam, ref_file, config, dbsnp,
                               region, out_file)
    if not file_exists(out_file):
        if not has_aligned_reads(align_bam, region):
            write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                params += ["-T", "UnifiedGenotyper",
                           "-o", tx_out_file,
                           "--genotype_likelihoods_model", "BOTH"]
                broad_runner.run_gatk(params)
    return out_file

def haplotype_caller(align_bam, ref_file, config, dbsnp=None,
                       region=None, out_file=None):
    """Call variation with GATK's HaplotypeCaller.

    This requires the full non open-source version of GATK.
    """
    broad_runner, params, out_file = \
        _shared_gatk_call_prep(align_bam, ref_file, config, dbsnp,
                               region, out_file)
    assert broad_runner.has_gatk_full(), \
        "Require full version of GATK 2.0 for haplotype based calling"
    if not file_exists(out_file):
        if not has_aligned_reads(align_bam, region):
            write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                params += ["-T", "HaplotypeCaller",
                           "-o", tx_out_file]
                broad_runner.run_gatk(params)
    return out_file

def write_empty_vcf(out_file):
    with open(out_file, "w") as out_handle:
        out_handle.write("##fileformat=VCFv4.1\n"
                         "## No variants; no reads aligned in region\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

# ## Utility functions for dealing with VCF files

def split_snps_indels(broad_runner, orig_file, ref_file):
    """Split a variant call file into SNPs and INDELs for processing.
    """
    base, ext = os.path.splitext(orig_file)
    snp_file = "{base}-snp{ext}".format(base=base, ext=ext)
    indel_file = "{base}-indel{ext}".format(base=base, ext=ext)
    params = ["-T", "SelectVariants",
              "-R", ref_file,
              "--variant", orig_file]
    for out_file, select_type in [(snp_file, ["SNP"]),
                                  (indel_file, ["INDEL", "MIXED", "MNP",
                                                "SYMBOLIC", "NO_VARIATION"])]:
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cur_params = params + ["--out", tx_out_file]
                for x in select_type:
                    cur_params += ["--selectTypeToInclude", x]
                broad_runner.run_gatk(cur_params)
    return snp_file, indel_file

def combine_variant_files(orig_files, out_file, ref_file, config,
                          quiet_out=True):
    """Combine multiple VCF files into a single output file.
    """
    broad_runner = broad.runner_from_config(config)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "CombineVariants",
                      "-R", ref_file,
                      "--out", tx_out_file]
            priority_order = []
            for orig_file in orig_files:
                name = os.path.splitext(os.path.basename(orig_file))[0]
                params.extend(["--variant:{name}".format(name=name), orig_file])
                priority_order.append(name)
            params.extend(["--rod_priority_list", ",".join(priority_order)])
            if quiet_out:
                params.extend(["--suppressCommandLineHeader", "--setKey", "null"])
            broad_runner.run_gatk(params)
    return out_file

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, config):
    """Filter variant calls using Variant Quality Score Recalibration.

    Newer GATK with Haplotype calling has combined SNP/indel filtering.
    """
    broad_runner = broad.runner_from_config(config)
    caller = config["algorithm"].get("variantcaller")
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    if caller in ["gatk-haplotype"] and cov_interval not in ["regional"]:
        return _variant_filtration_both(broad_runner, call_file, ref_file, vrn_files,
                                        config)
    # no additional filtration for callers that filter as part of call process
    elif caller in ["samtools"]:
        return call_file
    else:
        snp_file, indel_file = split_snps_indels(broad_runner, call_file, ref_file)
        snp_filter_file = _variant_filtration_snp(broad_runner, snp_file, ref_file,
                                                  vrn_files, config)
        indel_filter_file = _variant_filtration_indel(broad_runner, indel_file,
                                                      ref_file, vrn_files, config)
        orig_files = [snp_filter_file, indel_filter_file]
        out_file = "{base}combined.vcf".format(base=os.path.commonprefix(orig_files))
        return combine_variant_files(orig_files, out_file, ref_file, config)

def _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                         tranch_file, filter_type):
    """Apply recalibration details, returning filtered VCF file.
    """
    base, ext = os.path.splitext(snp_file)
    out_file = "{base}-{filter}filter{ext}".format(base=base, ext=ext,
                                                   filter=filter_type)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "ApplyRecalibration",
                      "-R", ref_file,
                      "--input", snp_file,
                      "--out", tx_out_file,
                      "--tranches_file", tranch_file,
                      "--recal_file", recal_file,
                      "--mode", filter_type]
            broad_runner.run_gatk(params)
    return out_file

def _shared_variant_filtration(filter_type, cov_interval, snp_file,
                               ref_file, vrn_files):
    """Share functionality for filtering variants.
    """
    recal_file = "{base}.recal".format(base = os.path.splitext(snp_file)[0])
    tranches_file = "{base}.tranches".format(base = os.path.splitext(snp_file)[0])
    params = ["-T", "VariantRecalibrator",
              "-R", ref_file,
              "--input", snp_file,
              "--mode", filter_type,
              "-an", "QD",
              "-an", "FS",
              "-an", "HaplotypeScore",
              "-an", "ReadPosRankSum"]
    if filter_type in ["SNP", "BOTH"]:
        params.extend(
            ["-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0",
             vrn_files.train_hapmap,
             "-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0",
             vrn_files.train_1000g_omni,
             "-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0",
             vrn_files.dbsnp,
              "-an", "MQRankSum",
              "-an", "MQ"])
    if filter_type in ["INDEL", "BOTH"]:
        assert vrn_files.train_indels, \
               "Need indel training file specified"
        params.extend(
            ["-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0",
             vrn_files.train_indels,])
    if cov_interval == "exome":
        params.extend(["--maxGaussians", "4", "--percentBadVariants", "0.05"])
    else:
        params.extend(["-an", "DP"])
    return params, recal_file, tranches_file

def variant_filtration_with_exp(broad_runner, snp_file, ref_file, filter_type,
                                expressions):
    """Perform hard filtering with GATK using JEXL expressions.

    Variant quality score recalibration will not work on some regions; it
    requires enough positions to train the model. This provides a general wrapper
    around GATK to do cutoff based filtering.
    """
    base, ext = os.path.splitext(snp_file)
    out_file = "{base}-filter{ftype}{ext}".format(base=base, ext=ext,
                                                  ftype=filter_type)
    if not file_exists(out_file):
        logger.info("Hard filtering %s with %s" % (snp_file, expressions))
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantFiltration",
                      "-R", ref_file,
                      "--out", tx_out_file,
                      "--variant", snp_file]
            for exp in expressions:
                params.extend(["--filterName", "GATKStandard{e}".format(e=exp.split()[0]),
                               "--filterExpression", exp])
            broad_runner.run_gatk(params)
    return out_file

# ## SNP specific variant filtration

def _variant_filtration_snp(broad_runner, snp_file, ref_file, vrn_files,
                            config):
    """Filter SNP variant calls using GATK best practice recommendations.
    """
    filter_type = "SNP"
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    params, recal_file, tranches_file = _shared_variant_filtration(
        filter_type, cov_interval, snp_file, ref_file, vrn_files)
    assert vrn_files.train_hapmap and vrn_files.train_1000g_omni, \
           "Need HapMap and 1000 genomes training files"
    filters = ["QD < 2.0", "MQ < 40.0", "FS > 60.0",
               "MQRankSum < -12.5", "ReadPosRankSum < -8.0"]
    # GATK Haplotype caller (v2.2) appears to have much larger HaplotypeScores
    # resulting in excessive filtering, so avoid this metric
    if variantcaller not in ["gatk-haplotype"]:
        filters.append("HaplotypeScore > 13.0")
    if cov_interval == "regional" or variantcaller == "freebayes":
        return variant_filtration_with_exp(broad_runner, snp_file, ref_file, filter_type,
                                           filters)
    else:
        # also check if we've failed recal and needed to do strict filtering
        filter_file = "{base}-filterSNP.vcf".format(base = os.path.splitext(snp_file)[0])
        if file_exists(filter_file):
            config["algorithm"]["coverage_interval"] = "regional"
            return _variant_filtration_snp(broad_runner, snp_file, ref_file, vrn_files,
                                           config)
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                try:
                    broad_runner.run_gatk(params)
                # Can fail to run if not enough values are present to train. Rerun with regional
                # filtration approach instead
                except:
                    config["algorithm"]["coverage_interval"] = "regional"
                    return _variant_filtration_snp(broad_runner, snp_file, ref_file, vrn_files,
                                                   config)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Indel specific variant filtration

def _variant_filtration_indel(broad_runner, snp_file, ref_file, vrn_files,
                              config):
    """Filter indel variant calls using GATK best practice recommendations.
    """
    filter_type = "INDEL"
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    params, recal_file, tranches_file = _shared_variant_filtration(
        filter_type, cov_interval, snp_file, ref_file, vrn_files)
    if cov_interval in ["exome", "regional"]:
        return variant_filtration_with_exp(broad_runner, snp_file, ref_file, filter_type,
                                           ["QD < 2.0", "ReadPosRankSum < -20.0", "FS > 200.0"])
    else:
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                broad_runner.run_gatk(params)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Variant filtration for combined indels and SNPs

def _variant_filtration_both(broad_runner, snp_file, ref_file, vrn_files,
                              config):
    """Filter SNP and indel variant calls using GATK best practice recommendations.
    """
    filter_type = "BOTH"
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    params, recal_file, tranches_file = _shared_variant_filtration(
        filter_type, cov_interval, snp_file, ref_file, vrn_files)
    if not file_exists(recal_file):
        with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
            params.extend(["--recal_file", tx_recal,
                           "--tranches_file", tx_tranches])
            broad_runner.run_gatk(params)
    return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                tranches_file, filter_type)

# ## Variant evaluation

def gatk_evaluate_variants(vcf_file, ref_file, config, dbsnp=None, intervals=None):
    """Evaluate variants, return SNP counts and Transition/Transversion ratios.
    """
    runner = broad.runner_from_config(config)
    eval_file = variant_eval(vcf_file, ref_file, dbsnp, intervals, runner)
    stats = _extract_eval_stats(eval_file)
    return _format_stats(stats['called'])

def _format_stats(stats):
    """Convert statistics into high level summary of major variables.
    """
    total = sum(itertools.chain.from_iterable(s.itervalues() for s in stats.itervalues()))
    if total > 0:
        dbsnp = sum(stats['known'].itervalues()) / float(total) * 100.0
    else:
        dbsnp = -1.0
    tv_dbsnp = stats['known']['tv']
    ti_dbsnp = stats['known']['ti']
    tv_novel = stats['novel']['tv']
    ti_novel = stats['novel']['ti']
    if tv_novel > 0 and tv_dbsnp > 0:
        titv_all = float(ti_novel + ti_dbsnp) / float(tv_novel + tv_dbsnp)
        titv_dbsnp = float(ti_dbsnp) / float(tv_dbsnp)
        titv_novel = float(ti_novel) / float(tv_novel)
    else:
        titv_all, titv_dbsnp, titv_novel = (-1.0, -1.0, -1.0)
    return dict(total=total, dbsnp_pct = dbsnp, titv_all=titv_all,
                titv_dbsnp=titv_dbsnp, titv_novel=titv_novel)

def _extract_eval_stats(eval_file):
    """Parse statistics of interest from GATK output file.
    """
    stats = dict()
    for snp_type in ['called', 'filtered']:
        stats[snp_type]  = dict()
        for dbsnp_type in ['known', 'novel']:
            stats[snp_type][dbsnp_type] = dict(ti=0, tv=0)
    for line in _eval_analysis_type(eval_file, "Ti/Tv Variant Evaluator"):
        if line[1:3] == ['dbsnp', 'eval']:
            snp_type = line[3]
            dbsnp_type = line[5]
            try:
                cur = stats[snp_type][dbsnp_type]
            except KeyError:
                cur = None
            if cur:
                stats[snp_type][dbsnp_type]["ti"] = int(line[6])
                stats[snp_type][dbsnp_type]["tv"] = int(line[7])
    return stats

def _eval_analysis_type(in_file, analysis_name):
    """Retrieve data lines associated with a particular analysis.
    """
    supported_versions = ["v0.2", "v1.0", "v1.1"]
    with open(in_file) as in_handle:
        # read until we reach the analysis
        for line in in_handle:
            if line.startswith(("##:GATKReport", "#:GATKReport")):
                version = line.split()[0].split(".", 1)[-1].split(":")[0]
                assert version in supported_versions, \
                       "Unexpected GATKReport version: {0}".format(version)
                if line.find(analysis_name) > 0:
                    break
        # read off header lines
        for _ in range(1):
            in_handle.next()
        # read the table until a blank line
        for line in in_handle:
            if not line.strip():
                break
            parts = line.rstrip("\n\r").split()
            yield parts

def variant_eval(vcf_in, ref_file, dbsnp, target_intervals, picard):
    """Evaluate variants in comparison with dbSNP reference.
    """
    out_file = "%s.eval" % os.path.splitext(vcf_in)[0]
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantEval",
                      "-R", ref_file,
                      "--eval", vcf_in,
                      "--dbsnp", dbsnp,
                      "-ST", "Filter",
                      "-o", tx_out_file,
                      "-l", "INFO",
                      "--doNotUseAllStandardModules",
                      "--evalModule", "CompOverlap",
                      "--evalModule", "CountVariants",
                      "--evalModule", "GenotypeConcordance",
                      "--evalModule", "TiTvVariantEvaluator",
                      "--evalModule", "ValidationReport",
                      "--stratificationModule", "Filter"]
            if target_intervals:
                # BED file target intervals are explicit with GATK 1.3
                # http://getsatisfaction.com/gsa/topics/
                # gatk_v1_3_and_bed_interval_file_must_be_parsed_through_tribble
                if _is_bed_file(target_intervals):
                    flag = "-L:bed"
                else:
                    flag = "-L"
                params.extend([flag, target_intervals])
            picard.run_gatk(params)
    return out_file

def _is_bed_file(fname):
    """Simple check if a file is in BED format.
    """
    if fname.lower().endswith(".bed"):
        return True
    with open(fname) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) > 3:
                    try:
                        int(parts[1])
                        int(parts[2])
                        return True
                    except ValueError:
                        pass
                break
    return False

# ## High level functionality to run genotyping in parallel

def _get_variantcaller(data):
    return data["config"]["algorithm"].get("variantcaller", "gatk")

def combine_multiple_callers(data):
    """Collapse together variant calls from multiple approaches into variants
    """
    by_bam = collections.defaultdict(list)
    for x in data:
        by_bam[x[0]["work_bam"]].append(x[0])
    out = []
    for grouped_calls in by_bam.itervalues():
        ready_calls = [{"variantcaller": _get_variantcaller(x),
                        "vrn_file": x["vrn_file"]}
                       for x in grouped_calls]
        final = grouped_calls[0]
        def orig_variantcaller_order(x):
            return final["config"]["algorithm"]["orig_variantcaller"].index(x["variantcaller"])
        if len(ready_calls) > 1:
            final["variants"] = sorted(ready_calls, key=orig_variantcaller_order)
        else:
            final["variants"] = ready_calls
        out.append([final])
    return out

def _handle_multiple_variantcallers(data):
    """Split samples that potentially require multiple variant calling approaches.
    """
    assert len(data) == 1
    callers = _get_variantcaller(data[0])
    if isinstance(callers, basestring):
        return [data]
    else:
        out = []
        for caller in callers:
            base = copy.deepcopy(data[0])
            base["config"]["algorithm"]["orig_variantcaller"] = \
              base["config"]["algorithm"]["variantcaller"]
            base["config"]["algorithm"]["variantcaller"] = caller
            out.append([base])
        return out

def parallel_variantcall(sample_info, parallel_fn):
    """Provide sample genotyping, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if x[0]["config"]["algorithm"]["snpcall"]:
            to_process.extend(_handle_multiple_variantcallers(x))
        else:
            finished.append(x)
    if len(to_process) > 0:
        split_fn = process_bam_by_chromosome("-variants.vcf", "work_bam",
                                             dir_ext_fn = _get_variantcaller)
        processed = parallel_split_combine(to_process, split_fn, parallel_fn,
                                           "variantcall_sample",
                                           "combine_variant_files",
                                           "vrn_file", ["sam_ref", "config"])
        finished.extend(processed)
    return finished

def variantcall_sample(data, region=None, out_file=None):
    """Parallel entry point for doing genotyping of a region of a sample.
    """
    from bcbio.variation import freebayes, cortex, samtools
    caller_fns = {"gatk": unified_genotyper,
                  "gatk-haplotype": haplotype_caller,
                  "freebayes": freebayes.run_freebayes,
                  "cortex": cortex.run_cortex,
                  "samtools": samtools.run_samtools}
    if data["config"]["algorithm"]["snpcall"]:
        sam_ref = data["sam_ref"]
        config = data["config"]
        caller_fn = caller_fns[config["algorithm"].get("variantcaller", "gatk")]
        data["vrn_file"] = caller_fn(data["work_bam"], sam_ref, config,
                                     configured_ref_file("dbsnp", config, sam_ref),
                                     region, out_file)
    return [data]


def split_vcf(in_file, config, out_dir=None):
    """
    split a VCF file into separate files by chromosome
    requires tabix to be installed

    """
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(in_file), "split")

    fasta_file = config["ref"]["fasta"]
    fasta_index = fasta_file + ".fai"
    samtools_path = config["program"].get("samtools", "samtools")
    tabix_path = config["program"].get("tabix", "tabix")

    if not file_exists(fasta_index):
        samtools = sh.Command(samtools_path)
        samtools.faidx(fasta_file)

    # if in_file is not compressed, compress it
    (_, ext) = os.path.splitext(in_file)
    if ext is not ".gz":
        gzip_file = in_file + ".gz"
        if not file_exists(gzip_file):
            sh.bgzip("-c", in_file, _out=gzip_file)
        in_file = gzip_file

    # create tabix index
    tabix_index(in_file)

    # find the chromosome names from the fasta index file
    chroms = str(sh.cut("-f1", fasta_index)).split()

    # make outfile from chromosome name
    def chr_out(chrom):
        out_file = replace_suffix(append_stem(in_file, chrom), ".vcf")
        return os.path.join(out_dir, os.path.basename(out_file))

    # run tabix to break up the vcf file
    def run_tabix(chrom):
        tabix = sh.Command(tabix_path)
        out_file = chr_out(chrom)
        if file_exists(out_file):
            return out_file
        with file_transaction(out_file) as tmp_out_file:
            tabix("-h", in_file, chrom, _out=tmp_out_file)
        return out_file

    out_files = map(run_tabix, chroms)
    return out_files


def tabix_index(in_file, preset="vcf", config=None):
    """
    index a file using tabix

    """
    if config:
        tabix_path = config["program"].get("tabix", "tabix")
    else:
        tabix_path = sh.which("tabix")
    tabix = sh.Command(tabix_path)
    out_file = in_file + ".tbi"
    if file_exists(out_file):
        return out_file
    tabix("-p", preset, in_file)

    return out_file
