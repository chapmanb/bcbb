"""
paver tasks related to project management

NOTE: these tasks rely on the options dictionary defined in paver. Make sure you have set the correct variables.
"""
import os
import contextlib
import glob

from paver.easy import *
from bcbio.paver.misc import process_list, process_dict

@contextlib.contextmanager
def cd_workdir(wd):
    orig_dir = os.getcwd()
    try:
        os.chdir(wd)
        yield
    finally:
        os.chdir(orig_dir)

# Utility functions
def setup(item):
    run_conf = os.path.join(options.dirs.data, item, "project_run_info.yaml")
    cl = ["setup_project_files.py", run_conf, os.path.join(options.dirs.data, item), '--project_dir', options.dirs.top]
    sh(" ".join(cl))

def _sbatch(cl, tmpl, **kw):
    kw.update(command_str = " ".join(cl))
    with open(os.path.join(options.dirs.sbatch, "%s.sh" % (kw['jobname'])), 'w') as out_handle:
        out_handle.write(tmpl.render(**kw))


def mako_to_rst(tmpl, **tmpl_kw):
    """Convert mako template tmpl to rst using template keywords **tmpl_kw"""
    outfile = tmpl.filename.replace(".mako", ".rst")
    with open(outfile, "w") as out_handle:
        out_handle.write(tmpl.render(**tmpl_kw))


def _exome_pipeline_cl(item):
    fc = glob.glob(os.path.join(options.dirs.intermediate, "*_*"))[0]
    wd = os.path.join(options.dirs.intermediate, fc)
    cl = ['exome_pipeline.py', os.path.join(options.dirs.git,"proj_conf.yaml"), 
          wd, os.path.join(options.dirs.data, item, "project_run_info.yaml"), 
          '--project_dir=%s' %(options.dirs.top) ]
    options.sbatch.update(jobname="%s_exomepipe" % (item), workdir=wd)
    _sbatch(cl, options.mako.sbatch, **options.sbatch)


# Tasks
@task
def setup_illumina_project():
    """Setup illumina project files"""
    options.log.info("Running bcbio.paver.project.setup_illumina_project")
    if len(options.illumina.flowcell_ids) > 0:
        process_list(setup, options.illumina.flowcell_ids)
    else:
        print >> sys.stderr, "No flowcell ids: skipping"


@task
def sbatch_exome_pipeline():
    """Make sbatch files for exome_pipeline.py
    See sbatch directory (options.dirs.sbatch) for output"""
    options.log.info("Running sbatch_exome_pipeline")
    with cd_workdir(options.dirs.top):
        process_list(_exome_pipeline_cl, options.illumina.flowcell_ids)

@task
def sbatch_exomeQC():
    """Make sbatch files for exomeQC.R
    """
    YAMLFILE = Template

# @task
# def make_exomeQC_sbatch(options):
#     """Make sbatch file for exomeQC.R
#     See %(sbatch) directory for output""" % (options.dirs.sbatch_dir)
#     with cd_workdir(options.dirs.sbatch_dir):
#         process_list(_make_exomeQC_cl, options.illumina.fcids)

