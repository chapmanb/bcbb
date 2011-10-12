"""Transfer raw files from finished NGS runs for backup and storage.
"""
import fabric.api as fabric

from bcbio.pipeline import log
from bcbio.log import create_log_handler
from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.transfer import remote_copy


def long_term_storage(remote_info, config_file):
    config = load_config(config_file)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        log.info("Copying run data over to remote storage: %s" %
        config["store_host"])
        log.debug("The contents from AMQP for this dataset are:\n %s" %
        remote_info)
        _copy_for_storage(remote_info, config)


def _copy_for_storage(remote_info, config):
    """Securely copy files from remote directory to the storage server.

    This requires ssh public keys to be setup so that no password entry
    is necessary, Fabric is used to manage setting up copies on the remote
    storage server.
    """
    base_dir = config["store_dir"]
    try:
        protocol = config["transfer_protocol"]
    except KeyError:
        protocol = None
        pass

    fabric.env.host_string = "%s@%s" % \
    (config["store_user"], config["store_host"])
    remote_copy(remote_info, base_dir, protocol)
