"""Testing using config files for transfer settings.
"""
import os
import sys
import subprocess
import fabric.api as fabric
import fabric.contrib.files as fabric_files
import time

def test_analyze_finished_sqn():
	"""Test running the script with the config files in the
	sub directory tests/test_transfer_data as input.

	NOTE: Requires running Galaxy to use, with correct values in the configs.
	"""
	config_dir = os.path.join(os.path.dirnam(__file__), "test_transfer_data")
	cl = ["analyze_finished_sqn.py", 
		  os.path.join(config_dir, "universe_wsgi.ini"),
		  os.path.join(config_dir, "post_process.yaml")]

	subprocess.check_call(cl)

def test_analyze_finished_sqn_transfer_info():
	"""Test running the script with the config files in the
	sub directory tests/test_transfer_data as input.

	NOTE: Requires running Galaxy to use, with correct values in the configs.
	"""
	config_dir = os.path.join(os.path.dirnam(__file__), "test_transfer_data")
	cl = ["analyze_finished_sqn.py", 
		  os.path.join(config_dir, "universe_wsgi.ini"),
		  os.path.join(config_dir, "post_process.yaml"),
		  os.path.join(config_dir, "transfer_info.yaml")]

	subprocess.check_call(cl)

# To be able to import functions from the scripts for testing.
sys.path.append(os.path.realpath("../scripts"))
sys.path.append(os.path.realpath(".."))

from analyze_finished_sqn import _remote_copy

fabric.env.key_filename = ["/Users/val/.ssh/local_ssh"]

def _remove_transferred_files(store_dir):
	"""Remove the files transferred in a previous test.
	"""
	copy_to = os.path.realpath("test_transfer_data/copy_to")
	with fabric.settings(host_string = "val@localhost"):
		fabric.run("rm -r %s/%s" % (copy_to, store_dir))

def test__remote_copy():
	"""Sets up dictionaries simulating loaded remote_info and config
	from various sources. Then test transferring files with the function
	using the standard setting.

	Note that there need to be passphrase free SSH between this machine
	and what is set up as remote_info["hostname"]. If remote_info["hostname]
	is "localhost"; SSH to localhost need to be set up to be passphrase free.
	"""
	store_dir = os.path.realpath("test_transfer_data/copy_to")

	config = {"analysis": {"store_dir": store_dir}}

	copy_dir = os.path.realpath("test_transfer_data/to_copy")

	remote_info = {}
	remote_info["directory"] = copy_dir
	remote_info["to_copy"] = ["file1", "file2", "file3"]
	remote_info["user"] = "val"
	remote_info["hostname"] = "localhost"

	# Generate test files
	test_data = {}
	for test_file in remote_info["to_copy"]:
		with open("%s/%s" % (copy_dir, test_file), "w") as file_to_write:
			test_data[test_file] = str(time.time())
			file_to_write.write(test_data[test_file])

	# Perform copy with settings
	with fabric.settings(host_string = "%s@%s" % (remote_info["user"], remote_info["hostname"])):
		
		if fabric_files.exists("%s/%s" % (store_dir, os.path.split(copy_dir)[1])):
			_remove_transferred_files(os.path.split(copy_dir)[1])

		# Copy
		_remote_copy(remote_info, config)

	# Check of the copy succeeded
	for test_file in remote_info["to_copy"]:
		with open("%s/%s/%s" % (store_dir, os.path.split(copy_dir)[1], test_file), "r") as copied_file:
			read_data = copied_file.read()
			assert read_data == test_data[test_file], "File copy failed"

def _test__remote_copy_scp():
	"""Sets up dictionaries simulating loaded remote_info and config
	from various sources. Then test transferring files with the function using
	the scp setting.

	Note that there need to be passphrase free SSH between this machine
	and what is set up as remote_info["hostname"]. If remote_info["hostname]
	is "localhost"; SSH to localhost need to be set up to be passphrase free.
	"""
	store_dir = os.path.realpath("test_transfer_data/copy_to")

	config = {"analysis": {"store_dir": store_dir}}

	copy_dir = os.path.realpath("test_transfer_data/to_copy")

	remote_info = {}
	remote_info["directory"] = copy_dir
	remote_info["to_copy"] = ["file1", "file2", "file3"]
	remote_info["user"] = "val"
	remote_info["hostname"] = "localhost"

	with fabric.settings(host_string = "%s@%s" % (remote_info["user"], remote_info["hostname"])):
		if fabric_files.exists("%s/%s" % (store_dir, os.path.split(copy_dir)[1])):
			_remove_transferred_files(os.path.split(copy_dir)[1])

		# Standard settings
		_remote_copy(remote_info, config, protocol = "scp")
