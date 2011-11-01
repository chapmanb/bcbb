"""Commandline interaction with SLURM schedulers.
"""
import re
import subprocess

_jobid_pat = re.compile("Submitted batch job (?P<jobid>\d+)")

def submit_job(scheduler_args, command):
    """Submit a job to the scheduler, returning the supplied job ID.
    """
    cl = ["sbatch"] + scheduler_args
    p = subprocess.Popen(cl,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (status,err) = p.communicate("#! /bin/sh\n%s" % command)
    p.close()
    match = _jobid_pat.search(status)
    return match.groups("jobid")[0]

def stop_job(jobid):
    cl = ["scancel", jobid]
    subprocess.check_call(cl)

def are_running(jobids):
    """Check if all of the submitted job IDs are running.
    """
    run_info = subprocess.check_output(["squeue"])
    running = []
    for parts in (l.split() for l in run_info.split("\n") if l.strip()):
        if len(parts) >= 5:
            pid, _, _, _, status = parts[:5]
            if status.upper() in ["R"]:
                running.append(pid)
    want_running = set(running).intersection(set(jobids))
    return len(want_running) == len(jobids)
