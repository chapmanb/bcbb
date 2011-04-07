"""Fabric deployment file to install genomic data on remote instances.

Designed to automatically download and manage biologically associated
data on cloud instances like Amazon EC2.

Fabric (http://docs.fabfile.org) is used to manage the automation of
a remote server.

Usage:
    fab -i key_file -H servername -f data_fabfile.py install_data
"""
import os
import operator
import socket
import glob
from contextlib import contextmanager

import yaml

from fabric.main import load_settings
from fabric.api import *
from fabric.contrib.files import *
try:
    import boto
except ImportError:
    boto = None

# -- Host specific setup for various groups of servers.

env.config_dir = os.path.join(os.path.dirname(__file__), "config")
env.remove_old_genomes = False

def setup_environment():
    """Setup environment with required data file locations.
    """
    _add_defaults()
    _amazon_ec2_environment()

def _amazon_ec2_environment():
    """Setup for a ubuntu amazon ec2 share; defaults if not set.

    Need to pass in host and private key file on commandline:
        -H hostname -i private_key_file
    """
    if not env.has_key("user") or not env.user:
        env.user = 'ubuntu'

def _add_defaults():
    """Defaults from fabricrc.txt file; loaded if not specified at commandline.
    """
    if not env.has_key("distribution"):
        config_file = os.path.join(env.config_dir, "fabricrc.txt")
        if os.path.exists(config_file):
            env.update(load_settings(config_file))

# -- Configuration for genomes to download and prepare

class _DownloadHelper:
    def __init__(self):
        self.config = {}

    def ucsc_name(self):
        return None

    def _exists(self, fname, seq_dir):
        """Check if a file exists in either download or final destination.
        """
        return exists(fname) or exists(os.path.join(seq_dir, fname))

class UCSCGenome(_DownloadHelper):
    def __init__(self, genome_name):
        _DownloadHelper.__init__(self)
        self._name = genome_name
        self._url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips" % \
                genome_name

    def ucsc_name(self):
        return self._name

    def download(self, seq_dir):
        zipped_file = None
        genome_file = "%s.fa" % self._name
        if not self._exists(genome_file, seq_dir):
            zipped_file = self._download_zip(seq_dir)
            if zipped_file.endswith(".tar.gz"):
                run("tar -xzpf %s" % zipped_file)
            elif zipped_file.endswith(".zip"):
                run("unzip %s" % zipped_file)
            elif zipped_file.endswith(".gz"):
                run("gunzip -c %s > out.fa" % zipped_file)
            else:
                raise ValueError("Do not know how to handle: %s" % zipped_file)
            tmp_file = genome_file.replace(".fa", ".txt")
            with settings(warn_only=True):
                result = run("ls *.fa")
            # some UCSC downloads have the files in multiple directories
            # mv them to the parent directory and delete the child directories
            #ignore_random = " -a \! -name '*_random.fa' -a \! -name 'chrUn*'" \
            #        "-a \! -name '*hap*.fa'"
            ignore_random = ""
            if result.failed:
                run("find . -name '*.fa'%s -exec mv {} . \;" % ignore_random)
                run("find . -type d -a \! -name '\.' | xargs rm -rf")
            result = run("find . -name '*.fa'%s" % ignore_random)
            result = result.split("\n")
            result.sort()
            run("cat %s > %s" % (" ".join(result), tmp_file))
            run("rm -f *.fa")
            run("mv %s %s" % (tmp_file, genome_file))
        return genome_file, [zipped_file]

    def _download_zip(self, seq_dir):
        for zipped_file in ["chromFa.tar.gz", "%s.fa.gz" % self._name,
                            "chromFa.zip"]:
            if not self._exists(zipped_file, seq_dir):
                with settings(warn_only=True):
                    result = run("wget %s/%s" % (self._url, zipped_file))
                if not result.failed:
                    break
            else:
                break
        return zipped_file

class NCBIRest(_DownloadHelper):
    """Retrieve files using the TogoWS REST server pointed at NCBI.
    """
    def __init__(self, name, refs):
        _DownloadHelper.__init__(self)
        self._name = name
        self._refs = refs
        self._base_url = "http://togows.dbcls.jp/entry/ncbi-nucleotide/%s.fasta"

    def download(self, seq_dir):
        genome_file = "%s.fa" % self._name
        if not self._exists(genome_file, seq_dir):
            for ref in self._refs:
                run("wget %s" % (self._base_url % ref))
                run("ls -l")
                run("sed -i.bak -r -e '/1/ s/^>.*$/>%s/g' %s.fasta" % (ref,
                    ref))
                # sed in Fabric does not cd properly?
                #sed('%s.fasta' % ref, '^>.*$', '>%s' % ref, '1')
            tmp_file = genome_file.replace(".fa", ".txt")
            run("cat *.fasta > %s" % tmp_file)
            run("rm -f *.fasta")
            run("rm -f *.bak")
            run("mv %s %s" % (tmp_file, genome_file))
        return genome_file, []

class EnsemblGenome(_DownloadHelper):
    """Retrieve genome FASTA files from Ensembl.

    ftp://ftp.ensemblgenomes.org/pub/plants/release-3/fasta/
    arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR9.55.dna.toplevel.fa.gz
    ftp://ftp.ensembl.org/pub/release-56/fasta/
    caenorhabditis_elegans/dna/Caenorhabditis_elegans.WS200.56.dna.toplevel.fa.gz
    """
    def __init__(self, ensembl_section, release_number, release2, organism,
            name, convert_to_ucsc=False):
        _DownloadHelper.__init__(self)
        if ensembl_section == "standard":
            url = "ftp://ftp.ensembl.org/pub/"
        else:
            url = "ftp://ftp.ensemblgenomes.org/pub/%s/" % ensembl_section
        url += "release-%s/fasta/%s/dna/" % (release_number, organism.lower())
        self._url = url
        release2 = ".%s" % release2 if release2 else ""
        self._get_file = "%s.%s%s.dna.toplevel.fa.gz" % (organism, name,
                release2)
        self._name = name
        self._convert_to_ucsc = convert_to_ucsc

    def download(self, seq_dir):
        genome_file = "%s.fa" % self._name
        if not self._exists(self._get_file, seq_dir):
            run("wget %s%s" % (self._url, self._get_file))
        if not self._exists(genome_file, seq_dir):
            run("gunzip -c %s > %s" % (self._get_file, genome_file))
        if self._convert_to_ucsc:
            #run("sed s/ / /g %s" % genome_file)
            raise NotImplementedError("Replace with chr")
        return genome_file, [self._get_file]

class BroadGenome(_DownloadHelper):
    """Retrieve genomes organized and sorted by Broad for use with GATK.
    """
    def __init__(self, target_fasta):
        _DownloadHelper.__init__(self)
        self._target = target_fasta
        self._ftp_url = "ftp://ftp.broadinstitute.org/pub/seq/references/"

    def download(self, seq_dir):
        if not self._exists(self._target, seq_dir):
            run("wget %s/%s" % (self._ftp_url, self._target))
        return self._target, []

GENOMES_SUPPORTED = [
           ("phiX174", "phix", NCBIRest("phix", ["NC_001422.1"])),
           ("Scerevisiae", "sacCer2", UCSCGenome("sacCer2")),
           ("Mmusculus", "mm9", UCSCGenome("mm9")),
           ("Mmusculus", "mm8", UCSCGenome("mm8")),
           ("Hsapiens", "hg18", UCSCGenome("hg18")),
           ("Hsapiens", "hg18-broad", BroadGenome("Homo_sapiens_assembly18.fasta")),
           ("Hsapiens", "GRCh37", BroadGenome("Homo_sapiens_assembly19.fasta")),
           ("Hsapiens", "hg19", UCSCGenome("hg19")),
           ("Rnorvegicus", "rn4", UCSCGenome("rn4")),
           ("Xtropicalis", "xenTro2", UCSCGenome("xenTro2")),
           ("Athaliana", "araTha_tair9", EnsemblGenome("plants", "6", "",
               "Arabidopsis_thaliana", "TAIR9")),
           ("Dmelanogaster", "dm3", UCSCGenome("dm3")),
           ("Celegans", "WS210", EnsemblGenome("standard", "60", "60",
               "Caenorhabditis_elegans", "WS210")),
           ("Mtuberculosis_H37Rv", "mycoTube_H37RV", NCBIRest("mycoTube_H37RV",
               ["NC_000962"])),
           ("Msmegmatis", "92", NCBIRest("92", ["NC_008596.1"])),
           ("Paeruginosa_UCBPP-PA14", "386", NCBIRest("386", ["CP000438.1"])),
           ("Ecoli", "eschColi_K12", NCBIRest("eschColi_K12", ["U00096.2"])),
           ("Amellifera_Honeybee", "apiMel3", UCSCGenome("apiMel3")),
           ("Cfamiliaris_Dog", "canFam2", UCSCGenome("canFam2")),
           ("Drerio_Zebrafish", "danRer6", UCSCGenome("danRer6")),
           ("Ecaballus_Horse", "equCab2", UCSCGenome("equCab2")),
           ("Fcatus_Cat", "felCat3", UCSCGenome("felCat3")),
           ("Ggallus_Chicken", "galGal3", UCSCGenome("galGal3")),
           ("Tguttata_Zebra_finch", "taeGut1", UCSCGenome("taeGut1")),
          ]

GENOME_INDEXES_SUPPORTED = ["bowtie", "bwa", "maq", "novoalign", "novoalign-cs",
                            "ucsc", "eland", "bfast", "arachne"]
DEFAULT_GENOME_INDEXES = ["seq"]
CONFIG_FILE = os.path.join(os.path.dirname(__file__), "config", "biodata.yaml")

# -- Fabric instructions

def install_data(config_file=CONFIG_FILE):
    """Main entry point for installing useful biological data.
    """
    _check_version()
    setup_environment()
    _data_uniref()
    genomes, genome_indexes = _get_genomes(config_file)
    _data_ngs_genomes(genomes, genome_indexes + DEFAULT_GENOME_INDEXES)
    lift_over_genomes = [g.ucsc_name() for (_, _, g) in genomes if g.ucsc_name()]
    _data_liftover(lift_over_genomes)

def install_data_s3(config_file=CONFIG_FILE):
    """Install data using pre-existing genomes present on Amazon s3.
    """
    _check_version()
    setup_environment()
    genomes, genome_indexes = _get_genomes(config_file)
    _download_genomes(genomes, genome_indexes + DEFAULT_GENOME_INDEXES)

def upload_s3(config_file=CONFIG_FILE):
    """Upload prepared genome files by identifier to Amazon s3 buckets.
    """
    if boto is None:
        raise ImportError("boto must be installed to upload to Amazon s3")
    if env.host != "localhost" and not env.host.startswith(socket.gethostname()):
        raise ValueError("Need to run S3 upload on a local machine")
    _check_version()
    setup_environment()
    genomes, genome_indexes = _get_genomes(config_file)
    _data_ngs_genomes(genomes, genome_indexes + DEFAULT_GENOME_INDEXES)
    _upload_genomes(genomes, genome_indexes + DEFAULT_GENOME_INDEXES)

def _check_version():
    version = env.version
    if int(version.split(".")[0]) < 1:
        raise NotImplementedError("Please install fabric version 1 or better")

def _get_genomes(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    genomes = []
    for g in config["genomes"]:
        ginfo = None
        for info in GENOMES_SUPPORTED:
            if info[1] == g["dbkey"]:
                ginfo = info
                break
        assert ginfo is not None, "Did not find download info for %s" % g["dbkey"]
        name, gid, manager = ginfo
        manager.config = g
        genomes.append((name, gid, manager))
    return genomes, config["genome_indexes"]

# == Decorators and context managers

def _if_installed(pname):
    """Run if the given program name is installed.
    """
    def argcatcher(func):
        def decorator(*args, **kwargs):
            with settings(
                    hide('warnings', 'running', 'stdout', 'stderr'),
                    warn_only=True):
                result = run(pname)
            if result.return_code not in [127]:
                return func(*args, **kwargs)
        return decorator
    return argcatcher

@contextmanager
def _make_tmp_dir():
    work_dir = os.path.join(env.data_files, "tmp")
    if not exists(work_dir):
        run("mkdir %s" % work_dir)
    yield work_dir
    if exists(work_dir):
        run("rm -rf %s" % work_dir)

# == NGS

def _data_ngs_genomes(genomes, genome_indexes):
    """Download and create index files for next generation genomes.
    """
    genome_dir = os.path.join(env.data_files, "genomes")
    if not exists(genome_dir):
        run('mkdir %s' % genome_dir)
    for organism, genome, manager in genomes:
        cur_dir = os.path.join(genome_dir, organism, genome)
        if not exists(cur_dir):
            run('mkdir -p %s' % cur_dir)
        with cd(cur_dir):
            if env.remove_old_genomes:
                _clean_genome_directory()
            seq_dir = 'seq'
            ref_file, base_zips = manager.download(seq_dir)
            ref_file = _move_seq_files(ref_file, base_zips, seq_dir)
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(cur_dir, ref_file, genome, cur_indexes, manager.config)

def _index_to_galaxy(work_dir, ref_file, gid, genome_indexes, config):
    """Index sequence files and update associated Galaxy loc files.
    """
    INDEX_FNS = {
        "seq" : _index_sam,
        "bwa" : _index_bwa,
        "bowtie": _index_bowtie,
        "maq": _index_maq,
        "novoalign": _index_novoalign,
        "novoalign_cs": _index_novoalign_cs,
        "ucsc": _index_twobit,
        "eland": _index_eland,
        "bfast": _index_bfast,
        "arachne": _index_arachne
        }
    indexes = {}
    with cd(work_dir):
        for idx in genome_indexes:
            indexes[idx] = INDEX_FNS[idx](ref_file)
    for ref_index_file, cur_index, prefix, new_style in [
          ("sam_fa_indices.loc", indexes.get("seq", None), "index", False),
          ("alignseq.loc", indexes.get("ucsc", None), "seq", False),
          ("twobit.loc", indexes.get("ucsc", None), "", False),
          ("bowtie_indices.loc", indexes.get("bowtie", None), "", True),
          ("bwa_index.loc", indexes.get("bwa", None), "", True)]:
        if cur_index:
            str_parts = _build_galaxy_loc_line(gid, os.path.join(work_dir, cur_index),
                                               config, prefix, new_style)
            _update_loc_file(ref_index_file, str_parts)

def _build_galaxy_loc_line(dbkey, file_path, config, prefix, new_style):
    """Prepare genome information to write to a Galaxy *.loc config file.
    """
    build_id = config.get("build_id", dbkey)
    display_name = config.get("name", dbkey)
    if new_style:
        str_parts = [build_id, dbkey, display_name, file_path]
    else:
        str_parts = [dbkey, file_path]
    if prefix:
        str_parts.insert(0, prefix)
    return str_parts

def _clean_genome_directory():
    """Remove any existing sequence information in the current directory.
    """
    for dirname in GENOME_INDEXES_SUPPORTED + DEFAULT_GENOME_INDEXES:
        if exists(dirname):
            run("rm -rf %s" % dirname)

def _move_seq_files(ref_file, base_zips, seq_dir):
    if not exists(seq_dir):
        run('mkdir %s' % seq_dir)
    for move_file in [ref_file] + base_zips:
        if exists(move_file):
            run("mv %s %s" % (move_file, seq_dir))
    path, fname = os.path.split(ref_file)
    moved_ref = os.path.join(path, seq_dir, fname)
    assert exists(moved_ref), moved_ref
    return moved_ref

def _update_loc_file(ref_file, line_parts):
    """Add a reference to the given genome to the base index file.
    """
    if env.galaxy_base is not None:
        tools_dir = os.path.join(env.galaxy_base, "tool-data")
        if not exists(tools_dir):
            conf_file = "tool_data_table_conf.xml"
            run("mkdir -p %s" % tools_dir)
            put(os.path.join("installed_files", conf_file),
                os.path.join(env.galaxy_base, conf_file))
        add_str = "\t".join(line_parts)
        with cd(tools_dir):
            if not exists(ref_file):
                run("touch %s" % ref_file)
            if not contains(ref_file, add_str):
                append(ref_file, add_str)

@_if_installed("faToTwoBit")
def _index_twobit(ref_file):
    """Index reference files using 2bit for random access.
    """
    dir_name = "ucsc"
    ref_base = os.path.splitext(os.path.split(ref_file)[-1])[0]
    out_file = "%s.2bit" % ref_base
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
    with cd(dir_name):
        if not exists(out_file):
            run("faToTwoBit %s %s" % (os.path.join(os.pardir, ref_file),
                out_file))
    return os.path.join(dir_name, out_file)

def _index_bowtie(ref_file):
    dir_name = "bowtie"
    ref_base = os.path.splitext(os.path.split(ref_file)[-1])[0]
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("bowtie-build -f %s %s" % (
                os.path.join(os.pardir, ref_file),
                ref_base))
    return os.path.join(dir_name, ref_base)

def _index_bwa(ref_file):
    dir_name = "bwa"
    local_ref = os.path.split(ref_file)[-1]
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("ln -s %s" % os.path.join(os.pardir, ref_file))
            with settings(warn_only=True):
                result = run("bwa index -a bwtsw %s" % local_ref)
            # work around a bug in bwa indexing for small files
            if result.failed:
                run("bwa index %s" % local_ref)
            run("rm -f %s" % local_ref)
    return os.path.join(dir_name, local_ref)

def _index_maq(ref_file):
    dir_name = "maq"
    local_ref = os.path.split(ref_file)[-1]
    binary_out = "%s.bfa" % os.path.splitext(local_ref)[0]
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("ln -s %s" % os.path.join(os.pardir, ref_file))
            run("maq fasta2bfa %s %s" % (local_ref,
                binary_out))
    return os.path.join(dir_name, binary_out)

@_if_installed("novoindex")
def _index_novoalign(ref_file):
    dir_name = "novoalign"
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    ref_file = os.path.join(os.pardir, ref_file)
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("novoindex %s %s" % (index_name, ref_file))

@_if_installed("novoalignCS")
def _index_novoalign_cs(ref_file):
    dir_name = "novoalign_cs"
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    ref_file = os.path.join(os.pardir, ref_file)
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("novoindex -c %s %s" % (index_name, ref_file))

def _index_sam(ref_file):
    (ref_dir, local_file) = os.path.split(ref_file)
    with cd(ref_dir):
        if not exists("%s.fai" % local_file):
            run("samtools faidx %s" % local_file)
    return ref_file

@_if_installed("MakeLookupTable")
def _index_arachne(ref_file):
    """Index for Broad's Arachne aligner.
    """
    dir_name = "arachne"
    ref_base = os.path.splitext(os.path.split(ref_file)[-1])[0]
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("ln -s %s" % os.path.join(os.pardir, ref_file))
            ref_file = os.path.split(ref_file)[-1]
            run("MakeLookupTable SOURCE=%s OUT_HEAD=%s" % (ref_file,
                ref_base))
            run("fastaHeaderSizes FASTA=%s HEADER_SIZES=%s.headerSizes" %
                    (ref_file, ref_file))
            #run("rm -f %s" % ref_file)
    return os.path.join(dir_name, ref_base)

@_if_installed("squashGenome")
def _index_eland(ref_file):
    """Index for Solexa's Eland aligner.

    This is nasty since Eland will choke on large files like the mm9 and h18
    genomes. It also has a restriction on only having 24 larger reference 
    files per directory. This indexes files with lots of shorter sequences (like
    xenopus) as one file, and splits up other files, removing random and other
    associated chromosomes to avoid going over the 24 file limit.
    """
    dir_name = "eland"
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        num_refs = run("grep '^>' %s | wc -l" % ref_file)
        # For a lot of reference sequences, Eland needs them in 1 file
        if int(num_refs) > 239:
            run("squashGenome %s %s" % (dir_name, ref_file))
        # For large reference sequences, squash fails and need them split up
        else:
            tmp_dir = "tmp_seqparts"
            run("mkdir %s" % tmp_dir)
            run("seqretsplit -sequence %s -osdirectory2 %s -outseq ." %
                    (ref_file, tmp_dir))
            with cd(tmp_dir):
                result = run("ls *.fasta")
                result = result.split("\n")
            seq_files = [os.path.join(tmp_dir, f) for f in result]
            run("squashGenome %s %s" % (dir_name, " ".join(seq_files)))
            run("rm -rf %s" % tmp_dir)
            # Eland can only handle up to 24 reference files in a directory
            # If we have more, remove any with *random* in the name to get
            # below. This sucks, but seemingly no way around it because
            # Eland will choke on large reference files
            if int(num_refs) > 24:
                with cd(dir_name):
                    for remove_re in ["*random*", "*_hap*", "chrun_*"]:
                        with settings(warn_only=True):
                            run("rm -f %s" % remove_re)
                    new_count = run("ls | wc -l")
                    # Human is still too big, need to remove chromosome M
                    if int(new_count) // 2 > 24:
                        with settings(warn_only=True):
                            run("rm -f chrm*")

# -- Genome upload and download to Amazon s3 buckets

def _download_genomes(genomes, genome_indexes):
    """Download a group of genomes from Amazon s3 bucket.
    """
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not exists(org_dir):
            run('mkdir -p %s' % org_dir)
        for idx in genome_indexes:
            with cd(org_dir):
                if not exists(idx):
                    url = "https://s3.amazonaws.com/biodata/genomes/%s-%s.tar.xz" % (gid, idx)
                    run("wget %s" % url)
                    run("tar -xJvpf %s" % os.path.basename(url))
                    run("rm -f %s" % os.path.basename(url))
        ref_file = os.path.join(org_dir, "seq", "%s.fa" % gid)
        assert exists(ref_file), ref_file
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(org_dir, ref_file, gid, cur_indexes, manager.config)

def _upload_genomes(genomes, genome_indexes):
    """Upload our configured genomes to Amazon s3 bucket.
    """
    conn = boto.connect_s3()
    bucket = conn.create_bucket("biodata")
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, _) in genomes:
        cur_dir = os.path.join(genome_dir, orgname, gid)
        _clean_directory(cur_dir, gid)
        for idx in genome_indexes:
            idx_dir = os.path.join(cur_dir, idx)
            tarball = _tar_directory(idx_dir, "%s-%s" % (gid, idx))
            _upload_to_s3(tarball, bucket)
    bucket.make_public()

def _upload_to_s3(tarball, bucket):
    """Upload the genome tarball to s3.
    """
    s3_key_name = os.path.join("genomes", os.path.basename(tarball))
    if not bucket.get_key(s3_key_name):
        gb_size = int(run("du -sm %s" % tarball).split()[0]) / 1000.0
        print "Uploading %s %.1fGb" % (s3_key_name, gb_size)
        if gb_size < 4.9:
            s3_key = bucket.new_key(s3_key_name)
            s3_key.set_contents_from_filename(tarball, reduced_redundancy=True)
            s3_key.set_acl("public-read")
        else:
            _large_file_upload(bucket, s3_key_name, tarball)
            s3_key = bucket.get_key(s3_key_name)
            s3_key.set_acl("public-read")

def _large_file_upload(bucket, s3_key_name, tarball):
    """Upload large files using Amazon's multipart upload functionality.
    """
    def split_file(in_file):
        prefix = os.path.join(os.path.dirname(in_file),
                              "S3PART%s" % (os.path.basename(s3_key_name).split(".")[0]))
        if not exists("%saa" % prefix):
            run("split -b250m %s %s" % (in_file, prefix))
        return sorted(glob.glob("%s*" % prefix))

    mp = bucket.initiate_multipart_upload(s3_key_name)
    for i, part in enumerate(split_file(tarball)):
        print "Transferring", part
        with open(part) as t_handle:
            mp.upload_part_from_file(t_handle, i+1)
        os.remove(part)
    mp.complete_upload()

def _tar_directory(dir, tar_name):
    """Create a tarball of the directory.
    """
    base_dir, tar_dir = os.path.split(dir)
    tarball = os.path.join(base_dir, "%s.tar.xz" % tar_name)
    if not exists(tarball):
        with cd(base_dir):
            run("tar -cJvpf %s %s" % (os.path.basename(tarball), tar_dir))
    return tarball

def _clean_directory(dir, gid):
    """Clean duplicate files from directories before tar and upload.
    """
    # get rid of softlinks
    bowtie_ln = os.path.join(dir, "bowtie", "%s.fa" % gid)
    maq_ln = os.path.join(dir, "maq", "%s.fa" % gid)
    for to_remove in [bowtie_ln, maq_ln]:
        if exists(to_remove):
            run("rm -f %s" % to_remove)
    # remove any downloaded original sequence files
    remove_exts = ["*.gz", "*.zip"]
    with cd(os.path.join(dir, "seq")):
        for rext in remove_exts:
            fnames = run("find . -name '%s'" % rext)
            for fname in (f.strip() for f in fnames.split("\n") if f.strip()):
                run("rm -f %s" % fname)

# == Liftover files

def _data_liftover(lift_over_genomes):
    """Download chain files for running liftOver.

    Does not install liftOver binaries automatically.
    """
    lo_dir = os.path.join(env.data_files, "liftOver")
    if not exists(lo_dir):
        run("mkdir %s" % lo_dir)
    lo_base_url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/liftOver/%s"
    lo_base_file = "%sTo%s.over.chain.gz"
    for g1 in lift_over_genomes:
        for g2 in [g for g in lift_over_genomes if g != g1]:
            g2u = g2[0].upper() + g2[1:]
            cur_file = lo_base_file % (g1, g2u)
            non_zip = os.path.splitext(cur_file)[0]
            worked = False
            with cd(lo_dir):
                if not exists(non_zip):
                    with settings(warn_only=True):
                        result = run("wget %s" % (lo_base_url % (g1, cur_file)))
                    # Lift over back and forths don't always exist
                    # Only move forward if we found the file
                    if not result.failed:
                        worked = True
                        run("gunzip %s" % cur_file)
            if worked:
                ref_parts = [g1, g2, os.path.join(lo_dir, non_zip)]
                _update_loc_file("liftOver.loc", ref_parts)

# == UniRef
def _data_uniref():
    """Retrieve and index UniRef databases for protein searches.

    http://www.ebi.ac.uk/uniref/

    These are currently indexed for FASTA searches. Are other indexes desired?
    Should this be separated out and organized by program like genome data?
    This should also check the release note and automatically download and
    replace older versions.
    """
    site = "ftp://ftp.uniprot.org"
    base_url = site + "/pub/databases/uniprot/" \
               "current_release/uniref/%s/%s"
    for uniref_db in ["uniref50", "uniref90", "uniref100"]:
        work_dir = os.path.join(env.data_files, "uniref", uniref_db)
        if not exists(work_dir):
            run("mkdir -p %s" % work_dir)
        base_work_url = base_url % (uniref_db, uniref_db)
        fasta_url = base_work_url + ".fasta.gz"
        base_file = os.path.splitext(os.path.basename(fasta_url))[0]
        with cd(work_dir):
            if not exists(base_file):
                run("wget -c %s" % fasta_url)
                run("gunzip %s" % os.path.basename(fasta_url))
                run("wget %s" % (base_work_url + ".release_note"))
        _index_blast_db(work_dir, base_file, "prot")

def _index_blast_db(work_dir, base_file, db_type):
    """Index a database using blast+ for similary searching.
    """
    type_to_ext = dict(prot = ("phr", "pal"), nucl = ("nhr", "nal"))
    db_name = os.path.splitext(base_file)[0]
    with cd(work_dir):
        if not reduce(operator.or_,
            (exists("%s.%s" % (db_name, ext)) for ext in type_to_ext[db_type])):
            run("makeblastdb -in %s -dbtype %s -out %s" %
                    (base_file, db_type, db_name))


# == Not used -- takes up too much space and time to index

def _index_bfast(ref_file):
    """Indexes bfast in color and nucleotide space for longer reads.

    This preps for 40+bp sized reads, which is bfast's strength.
    """
    dir_name = "bfast"
    window_size = 14
    bfast_nt_masks = [
   "1111111111111111111111",
   "1111101110111010100101011011111",
   "1011110101101001011000011010001111111",
   "10111001101001100100111101010001011111",
   "11111011011101111011111111",
   "111111100101001000101111101110111",
   "11110101110010100010101101010111111",
   "111101101011011001100000101101001011101",
   "1111011010001000110101100101100110100111",
   "1111010010110110101110010110111011",
    ]
    bfast_color_masks = [
    "1111111111111111111111",
    "111110100111110011111111111",
    "10111111011001100011111000111111",
    "1111111100101111000001100011111011",
    "111111110001111110011111111",
    "11111011010011000011000110011111111",
    "1111111111110011101111111",
    "111011000011111111001111011111",
    "1110110001011010011100101111101111",
    "111111001000110001011100110001100011111",
    ]
    local_ref = os.path.split(ref_file)[-1]
    if not exists(dir_name):
        run("mkdir %s" % dir_name)
        with cd(dir_name):
            run("ln -s %s" % os.path.join(os.pardir, ref_file))
            # nucleotide space
            run("bfast fasta2brg -f %s -A 0" % local_ref)
            for i, mask in enumerate(bfast_nt_masks):
                run("bfast index -d 1 -n 4 -f %s -A 0 -m %s -w %s -i %s" %
                        (local_ref, mask, window_size, i + 1))
            # colorspace
            run("bfast fasta2brg -f %s -A 1" % local_ref)
            for i, mask in enumerate(bfast_color_masks):
                run("bfast index -d 1 -n 4 -f %s -A 1 -m %s -w %s -i %s" %
                        (local_ref, mask, window_size, i + 1))
