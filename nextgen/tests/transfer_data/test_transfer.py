"""Testing using config files for transfer settings.
"""
import os
import fabric.api as fabric
import fabric.contrib.files as fabric_files
import time
from bcbio.pipeline import log
from bcbio.pipeline.config_loader import load_config

from bcbio.pipeline.toplevel import _copy_from_sequencer
from bcbio.pipeline.storage import _copy_for_storage


def _remove_transferred_files(remote_info, config):
    """Remove the files transferred in a previous test.
    """
    copy_to = os.path.realpath("../transfer_data/copy_to")
    with fabric.settings(host_string="%s@%s" % \
         (config["store_user"], config["store_host"])):
        rm_str = "rm -r %s/%s" % \
         (copy_to, os.path.split(remote_info["directory"])[1])
        log.debug(rm_str)
        fabric.run(rm_str)


def perform_transfer(transfer_function, protocol_config, \
        remove_before_copy=True, should_overwrite=False):
    """Sets up dictionaries simulating loaded remote_info and config
    from various sources. Then test transferring files with the function
    using the standard setting.

    Note that there need to be passphrase free SSH between this machine
    and what is set up as store_host in post_process.yaml.

    The file post_process.yaml need the information "store_user" and
    "store_host".

    transfer_function : function - The function to use for transferring
    the files.
        This function should accept the two dictionaries config and
        remote_info as parameters.
    protocol_config : dictionary - containing the key "transfer_protocol"
    which give the string
    of a transfer protocol.
    """
    config = load_config("../data/automated/post_process.yaml")

    store_dir = os.path.realpath("../transfer_data/copy_to")
    config["store_dir"] = store_dir

    config.update(protocol_config)

    copy_dir = os.path.realpath("../transfer_data/to_copy")

    remote_info = {}
    remote_info["directory"] = copy_dir
    remote_info["to_copy"] = ["file1", "file2", "file3", "dir1"]
    remote_info["user"] = config["store_user"]
    remote_info["hostname"] = config["store_host"]

    # Generate test files
    for file_dir in store_dir, copy_dir:
        if not os.path.isdir(file_dir):
            os.mkdir(file_dir)

    test_data = {}
    for test_file in remote_info["to_copy"][:3]:
        test_file_path = "%s/%s" % (copy_dir, test_file)
        if not os.path.isdir(test_file_path):
            with open(test_file_path, 'w+') as file_to_write:
                # We just use the current processor time as test data, the
                # important part is that it will be different enough between
                # test just so we know we are not comparing with files copied
                # in a previous test during the assertion.
                test_data[test_file] = str(time.clock())
                file_to_write.write(test_data[test_file])

    if not os.path.isdir("%s/%s" % (copy_dir, remote_info["to_copy"][3])):
        os.mkdir("%s/%s" % (copy_dir, remote_info["to_copy"][3]))

    # Perform copy with settings
    with fabric.settings(host_string="%s@%s" % \
         (remote_info["user"], remote_info["hostname"])):

        if fabric_files.exists("%s/%s" % \
           (store_dir, os.path.split(copy_dir)[1])) and remove_before_copy:
            _remove_transferred_files(remote_info, config)

        # Copy
        transfer_function(remote_info, config)

    # Check of the copy succeeded
    for test_file in remote_info["to_copy"][:3]:
        test_file_path = "%s/%s/%s" % \
                         (store_dir, os.path.split(copy_dir)[1], test_file)
        # Did the files get copied correctly
        assert os.path.isfile(test_file_path), "File not copied: %s" % test_file
        if os.path.isfile(test_file_path):
            with open(test_file_path, "r") as copied_file:
                read_data = copied_file.read()
                fail_string = "File copy failed: %s " % read_data + \
                "doesn't equal %s (for %s). " % \
                (test_data[test_file], test_file_path) + \
                "Remove is %s and Overwrite is %s." % \
                (str(remove_before_copy), str(should_overwrite))
                # Assertion that passes when:
                #  - The files got copied if we removed the old files in the
                # target directory before the copy.
                #  - The new files did not replace the old files in the
                # target directory if they where not supposed to overwrite.
                #  - The new files did replace the old files in the target
                # directory of we specified that this should happen.
                assert (read_data == test_data[test_file]) == \
                (remove_before_copy or should_overwrite), fail_string
        # Did the directories get copied correcty
        if os.path.isdir(test_file_path):
            # Not tested for yet
            pass


def test__copy_for_storage():
    """Test using the copy function without any specification
    as to how to do it.
    """
    config = {}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False)


def test__copy_for_storage_scp():
    """Test using the copy function with scp.
    """
    config = {"transfer_protocol": "scp"}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False)


def test__copy_for_storage_rsync():
    """Test using the copy function with rsync.
    """
    config = {"transfer_protocol": "rsync"}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False, should_overwrite=True)

# NOTE: rdiff-backup not available
# def test__copy_for_storage_rdiff_backup():
#     """Test using the copy function with rdiff-backup.
#     """
#     config = {"transfer_protocol": "rdiff-backup"}
#     perform_transfer(_copy_for_storage, config)
#     perform_transfer(_copy_for_storage, config,
#     remove_before_copy=False, should_overwrite=True)


def test__copy_from_sequencer():
    """Test using the copy function without any specification
    as to how to do it.
    """
    config = {}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False)


def test__copy_from_sequencer_scp():
    """Test using the copy function with scp.
    """
    config = {"transfer_protocol": "scp"}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False)


def test__copy_from_sequencer_rsync():
    """Test using the copy function with rsync.
    """
    config = {"transfer_protocol": "rsync"}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False, should_overwrite=True)

# NOTE: rdiff-backup not available
# def test__copy_from_sequencer_rdiff_backup():
#     """Test using the copy function with rdiff-backup.
#     """
#     config = {"transfer_protocol": "rdiff-backup"}
#     perform_transfer(_copy_from_sequencer, config)
#     perform_transfer(_copy_from_sequencer, config,
#     remove_before_copy=False, should_overwrite=True)
