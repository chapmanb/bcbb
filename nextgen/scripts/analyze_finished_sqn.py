"""Server which listens for finished NGS runs, processing for upload to galaxy.

Usage:
    analyze_finished_sqn.py <Galaxy config file> <Post-processing config file>

The server can run the copy and processing on a remote host by setting the
analysis user and host parameters in your post_process.yaml file.
ssh keys need to be configured to allow passwordless login between
the machine running this server and the analysis host.

Need to configure the RabbitMQ server with:

    rabbitmqctl add_user galaxy password
    rabbitmqctl add_vhost galaxy_messaging_engine
    rabbitmqctl set_permissions -p galaxy_messaging_engine galaxy '.*' '.*' '.*'
"""
import os
import re
import sys
import ConfigParser
import json
import subprocess
import contextlib
import logbook

import yaml
from amqplib import client_0_8 as amqp
import fabric.api as fabric
import fabric.contrib.files as fabric_files

from bcbio.log import create_log_handler

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(galaxy_config, processing_config):
    amqp_config = _read_amqp_config(galaxy_config)
    with open(processing_config) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, LOG_NAME)
    process_tag = config["msg_process_tag"]
    handlers = [(process_tag,
                 analysis_handler(config, process_tag, processing_config))]
    with log_handler.applicationbound():
        message_reader(handlers, amqp_config)

def copy_and_analyze(remote_info, config, config_file):
    """Remote copy an output directory, process it, and upload to Galaxy.
    """
    log.debug("Remote host information: %s" % remote_info)
    a_host_str, a_shell, c_host_str = _config_hosts(config)

    with fabric.settings(host_string=c_host_str):
        fc_dir = _remote_copy(remote_info, config)

    with fabric.settings(host_string=a_host_str, shell=a_shell):
        _analyze_and_upload(remote_info, config, config_file, fc_dir)

def _analyze_and_upload(remote_info, config, config_file, fc_dir):
    if config["analysis"].get("config_file", None):
        config_file = config["analysis"]["config_file"]
    elif not config_file.startswith("/"):
        config_file = os.path.join(os.getcwd(), config_file)

    # Converted from an Illumina/Genesifter SampleSheet.csv
    run_yaml = os.path.join(config["analysis"]["store_dir"],
                                os.path.basename(fc_dir), "run_info.yaml")
    
    # TODO: Find a solution for PUSH architecture (check if *remote* file exists too)
    # fabric_files.exists(run_yaml)
    #if not os.path.exists(run_yaml):
    #    run_yaml = None

    analysis_dir = os.path.join(config["analysis"]["base_dir"],
                                os.path.basename(remote_info["directory"]))
    if not fabric_files.exists(analysis_dir):
        fabric.run("mkdir %s" % analysis_dir)
    with fabric.cd(analysis_dir):
        cl = [config["analysis"]["process_program"], config_file, fc_dir]
        if run_yaml:
            cl.append(run_yaml)
        fabric.run(" ".join(cl))
    cl = [config["analysis"]["upload_program"], config_file, fc_dir, analysis_dir]
    fabric.run(" ".join(cl))

def _config_hosts(config):
    """Retrieve configured machines to perform analysis and copy on.
    """
    user = config["analysis"].get("user", None)
    host = config["analysis"].get("host", None)
    shell = config["analysis"].get("login_shell", None)
    if not user or not host:
        user = os.environ["USER"]
        host = re.sub(r'\..*', '', os.uname()[1])
    if not shell:
        shell = os.environ["SHELL"]
    copy_user = config["analysis"].get("copy_user", None)
    copy_host = config["analysis"].get("copy_host", None)
    if not copy_user or not copy_host:
        copy_user, copy_host = (user, host)
    analysis_host_str = "%s@%s" % (user, host)
    analysis_shell = "%s -i -l -c" % shell
    copy_host_str = "%s@%s" % (copy_user, copy_host)
    return analysis_host_str, analysis_shell, copy_host_str

def _remote_copy(remote_info, config):
    """Securely copy files from remote directory to the processing server.

    This requires ssh public keys to be setup so that no password entry
    is necessary.
    """
    fc_dir = os.path.join(config["analysis"]["store_dir"],
                          os.path.basename(remote_info['directory']))
    log.info("Copying analysis files to %s" % fc_dir)
    if not fabric_files.exists(fc_dir):
        fabric.run("mkdir %s" % fc_dir)
    for fcopy in remote_info['to_copy']:
        target_loc = os.path.join(fc_dir, fcopy)
        if not fabric_files.exists(target_loc):
            target_dir = os.path.dirname(target_loc)
            if not fabric_files.exists(target_dir):
                fabric.run("mkdir -p %s" % target_dir)
            cl = ["scp", "-r", "%s@%s:%s/%s" %
                  (remote_info["user"], remote_info["hostname"],
                   remote_info["directory"], fcopy),
                  target_loc]
            fabric.run(" ".join(cl))
    log.info("Analysis files copied")
    return fc_dir

def analysis_handler(processing_config, tag_name, config_file):
    def receive_msg(msg):
        if msg.properties['application_headers'].get('msg_type') == tag_name:
            copy_and_analyze(json.loads(msg.body), processing_config,
                config_file)
    return receive_msg

def message_reader(handlers, config):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    conn = amqp.Connection(host=config['host'] + ":" + config['port'],
                           userid=config['userid'], password=config['password'],
                           virtual_host=config['virtual_host'], insist=False)
    chan = conn.channel()
    for tag_name, handler in handlers:

        chan.queue_declare(queue=tag_name, exclusive=False, auto_delete=False,
                durable=True)
        try:
            chan.exchange_declare(exchange=config['exchange'], type="fanout", durable=True,
                    auto_delete=False)
        except amqp.exceptions.AMQPChannelException:
            chan.exchange_delete(exchange=config['exchange'])
            chan.exchange_declare(exchange=config['exchange'], type="fanout", durable=True,
                    auto_delete=False)
        chan.queue_bind(queue=tag_name, exchange=config['exchange'],
                        routing_key=config['routing_key'])
        log.debug("AMQP: Consuming %s" % tag_name)
        chan.basic_consume(queue=tag_name, no_ack=True,
                           callback=handler, consumer_tag=tag_name)

    log.debug("AMQP: Waiting to consume message")
    while True:
        chan.wait()
    for (tag_name, _) in handlers:
        chan.basic_cancel(tag_name)
    chan.close()
    conn.close()

def _read_amqp_config(galaxy_config):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    config = ConfigParser.ConfigParser()
    config.read(galaxy_config)
    amqp_config = {}
    for option in config.options("galaxy_amqp"):
        amqp_config[option] = config.get("galaxy_amqp", option)
    return amqp_config

if __name__ == "__main__":
    main(*sys.argv[1:])
