"""Testing using config files for transfer settings.
"""
import os
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

def test__remote_copy():
	"""Sets up dictionaries simulating loaded remote_info and config
	from various sources. Then test transferring files with the function.
	"""
	config = {'analysis': {'store_dir': '/Users/val/pipeline_test/store_dir'}}
	remote_info = {}
	remote_info["directory"] = "/Users/val/Documents/bcbb/nextgen/tests/test_transfer_data/to_copy"
	remote_info["to_copy"] = ["file1", "file2", "file3"]
	remote_info["user"] = "val"
	remote_info["hostname"] = "localhost"

	_remote_copy(remote_info, config)