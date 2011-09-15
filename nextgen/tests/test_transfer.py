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

def get_transfer_function(transfer_config):
	"""Returns a function to use for transfer where we will have set the
	parameter protocol=setting

	transfer_config : dictionary - The dictionary containing transfer
		configurations. For example which protocol to use.
	"""
	def transfer_function(remote_info, config):
		_remote_copy(remote_info, config, transfer_config = transfer_config)

	return transfer_function

def perform__remote_copy_test(transfer_function):
	"""Sets up dictionaries simulating loaded remote_info and config
	from various sources. Then test transferring files with the function
	using the standard setting.

	Note that there need to be passphrase free SSH between this machine
	and what is set up as remote_info["hostname"]. If remote_info["hostname]
	is "localhost"; SSH to localhost need to be set up to be passphrase free.

	transfer_function : function - The function to use for transferring
	the files.
		This function should accept the two dictionaries config and
		remote_info as parameters.
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
			# We just the current processor time as test data, the important
			# part is that it will be different enough between tests just
			# so we know we are not comparing with files copied in a 
			# previous test during the assertion.
			test_data[test_file] = str(time.clock())
			file_to_write.write(test_data[test_file])

	# Perform copy with settings
	with fabric.settings(host_string = "%s@%s" % (remote_info["user"], remote_info["hostname"])):
		
		if fabric_files.exists("%s/%s" % (store_dir, os.path.split(copy_dir)[1])):
			_remove_transferred_files(os.path.split(copy_dir)[1])

		# Copy
		transfer_function(remote_info, config)

	# Check of the copy succeeded
	for test_file in remote_info["to_copy"]:
		with open("%s/%s/%s" % (store_dir, os.path.split(copy_dir)[1], test_file), "r") as copied_file:
			read_data = copied_file.read()
			assert read_data == test_data[test_file], "File copy failed"

def test__remote_copy():
	"""Test using the copy function without any specification
	as to how to do it.
	"""
	perform__remote_copy_test(_remote_copy)

def test__remote_copy_scp():
	"""Test using the copy function with scp.
	"""
	transfer_config = {"transfer_protocol" : "scp"}
	copy_function = get_transfer_function(transfer_config)
	perform__remote_copy_test(copy_function)

def test__remote_copy_rsync():
	"""Test using the copy function with rsync.
	"""
	transfer_config = {"transfer_protocol" : "rsync"}
	copy_function = get_transfer_function(transfer_config)
	perform__remote_copy_test(copy_function)

def test__remote_copy_rdiff_backup():
	"""Test using the copy function with rdiff-backup.
	"""
	transfer_config = {"transfer_protocol" : "rdiff-backup"}
	copy_function = get_transfer_function(transfer_config)
	perform__remote_copy_test(copy_function)
