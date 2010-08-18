"""Install next gen sequencing analysis tools not currently packaged.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir

@_if_not_installed("faToTwoBit")
def install_ucsc_tools(env):
    """Install useful executables from UCSC.

    todo: install from source to handle 32bit and get more programs
    http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
    """
    tools = ["liftOver", "faToTwoBit", "bedToBigBed",
             "bigBedInfo", "bigBedSummary", "bigBedToBed",
             "bigWigInfo", "bigWigSummary", "bigWigToBedGraph", "bigWigToWig",
             "fetchChromSizes", "wigToBigWig"]
    url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
    install_dir = os.path.join(env.system_install, "bin")
    for tool in tools:
        with cd(install_dir):
            if not exists(tool):
                sudo("wget %s%s" % (url, tool))
                sudo("chmod a+rwx %s" % tool)

@_if_not_installed("bowtie")
def install_bowtie(env):
    """Install the bowtie short read aligner.
    """
    version = "0.12.5"
    mirror_info = "?use_mirror=cdnetworks-us-1"
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie/%s/" \
          "bowtie-%s-src.zip" % (version, version)
    install_dir = os.path.join(env.system_install, "bin")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s" % (url, mirror_info))
            run("unzip %s" % os.path.split(url)[-1])
            with cd("bowtie-%s" % version):
                run("make")
                for fname in run("find -perm -100 -name 'bowtie*'").split("\n"):
                    sudo("mv -f %s %s" % (fname, install_dir))

@_if_not_installed("bwa")
def install_bwa(env):
    version = "0.5.7"
    mirror_info = "?use_mirror=cdnetworks-us-1"
    url = "http://downloads.sourceforge.net/project/bio-bwa/bwa-%s.tar.bz2" % (
            version)
    install_dir = os.path.join(env.system_install, "bin")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s" % (url, mirror_info))
            run("tar -xjvpf %s" % (os.path.split(url)[-1]))
            with cd("bwa-%s" % version):
                arch = run("uname -m")
                # if not 64bit, remove the appropriate flag
                if arch.find("x86_64") == -1:
                    run("sed -i.bak -r -e 's/-O2 -m64/-O2/g' Makefile")
                run("make")
                sudo("mv bwa %s" % install_dir)
                sudo("mv solid2fastq.pl %s" % install_dir)
                sudo("mv qualfa2fq.pl %s" % install_dir)

@_if_not_installed("fastq_quality_boxplot_graph.sh")
def install_fastx_toolkit(env):
    version = "0.0.13"
    gtext_version = "0.6"
    url_base = "http://hannonlab.cshl.edu/fastx_toolkit/"
    fastx_url = "%sfastx_toolkit-%s.tar.bz2" % (url_base, version)
    gtext_url = "%slibgtextutils-%s.tar.bz2" % (url_base, gtext_version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % gtext_url)
            run("tar -xjvpf %s" % (os.path.split(gtext_url)[-1]))
            with cd("libgtextutils-%s" % gtext_version):
                run("./configure --prefix=%s" % (env.system_install))
                run("make")
                sudo("make install")
            run("wget %s" % fastx_url)
            run("tar -xjvpf %s" % os.path.split(fastx_url)[-1])
            with cd("fastx_toolkit-%s" % version):
                run("./configure --prefix=%s" % (env.system_install))
                run("make")
                sudo("make install")

@_if_not_installed("bfast")
def install_bfast(env):
    version = "0.6.4"
    vext = "d"
    url = "http://downloads.sourceforge.net/project/bfast/bfast/%s/bfast-%s%s.tar.gz"\
            % (version, version, vext)
    install_dir = os.path.join(env.system_install, "bin")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % (url))
            run("tar -xzvpf %s" % (os.path.split(url)[-1]))
            with cd("bfast-%s%s" % (version, vext)):
                run("./configure --prefix=%s" % (install_dir))
                run("make")
                sudo("make install")

def _symlinked_java_version_dir(pname, version):
    base_dir = os.path.join(env.system_install, "share", "java", pname)
    install_dir = "%s-%s" % (install_dir, version)
    if not exists(install_dir):
        sudo("mkdir -p %s" % install_dir)
        if exists(base_dir):
            sudo("rm -f %s" % base_dir)
        sudo("ln -s %s %s" % (install_dir, base_dir))
        return install_dir
    return None

def install_picard(env):
    version = "1.26"
    url = "http://downloads.sourceforge.net/project/picard/" \
          "picard-tools/%s/picard-tools-%s.zip" % (version, version)
    install_dir = _symlinked_java_version_dir("picard", version)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget %s" % (url))
                run("unzip %s" % os.path.basename(url))
                with cd(os.path.splitext(os.path.basename(url))[0]):
                    sudo("mv *.jar %s" % install_dir)

def install_gatk(env):
    version = "1.0.3857"
    ext = ".tar.bz2"
    url = "ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/"\
          "GenomeAnalysisTK-%s%s" % (version, ext)
    install_dir = _symlinked_java_version_dir("gatk", version)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget %s" % (url))
                run("tar -xjvpf %s" % os.path.basename(url))
                with cd(os.path.basename(url).replace(ext, "")):
                    sudo("mv *.jar %s" % install_dir)
