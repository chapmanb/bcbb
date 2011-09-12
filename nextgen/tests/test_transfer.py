"""Testing using config files for transfer settings.
"""
import os
import sys
import subprocess

def test_analyze_finished_sqn():
	"""Test running the script with the config files in the
	sub directory tests/test_transfer_data as input.

	Requires running Galaxy to use, with correct values in the configs.
	"""
	config_dir = os.path.join(os.path.dirnam(__file__), "test_transfer_data")
	cl = ["analyze_finished_sqn.py", 
		  os.path.join(config_dir, "universe_wsgi.ini"),
		  os.path.join(config_dir, "post_process.yaml")]

	subprocess.check_call(cl)

def test_analyze_finished_sqn_transfer_info():
	"""Test running the script with the config files in the
	sub directory tests/test_transfer_data as input.

	Requires running Galaxy to use, with correct values in the configs.
	"""
	config_dir = os.path.join(os.path.dirnam(__file__), "test_transfer_data")
	cl = ["analyze_finished_sqn.py", 
		  os.path.join(config_dir, "universe_wsgi.ini"),
		  os.path.join(config_dir, "post_process.yaml"),
		  os.path.join(config_dir, "transfer_info.yaml")]

	subprocess.check_call(cl)

# To be able to import functions from the scripts for testing.
sys.path.append("/Users/val/Documents/bcbb/nextgen/scripts")
sys.path.append("/Users/val/Documents/bcbb/nextgen")

from analyze_finished_sqn import _remote_copy
import fabric.api as fabric

def test__remote_copy():
	"""Sets up dictionaries simulating loaded remote_info and config
	from various sources. Then test transferring files with the function.

	Note that there need to be passphrase free SSH between this machine
	and what is set up as remote_info["hostname"]. If remote_info["hostname]
	is "localhost"; SSH to localhost need to be set up to be passphrase free.
	"""
	config = {"analysis": {"store_dir": "/Users/val/pipeline_test/store_dir"}}
	remote_info = {}
	remote_info["directory"] = "/Users/val/Documents/bcbb/nextgen/tests/test_transfer_data/to_copy"
	remote_info["to_copy"] = ["file1", "file2", "file3"]
	remote_info["user"] = "val"
	remote_info["hostname"] = "localhost"

	with fabric.settings(host_string = "%s@%s" % (remote_info["user"], remote_info["hostname"])):
		_remote_copy(remote_info, config)
