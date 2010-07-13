"""Server which listens for finished NGS runs, processing for upload to galaxy.

Usage:
    analyze_finished_sqn.py <Galaxy config file> <Post-processing config file>

Need to configure the RabbitMQ server with:

    rabbitmqctl add_user galaxy password
    rabbitmqctl add_vhost galaxy_messaging_engine
    rabbitmqctl set_permissions -p galaxy_messaging_engine galaxy '.*' '.*' '.*'
"""
import sys
import ConfigParser
import json
import subprocess

import yaml
from amqplib import client_0_8 as amqp

def main(galaxy_config, processing_config):
    amqp_config = _read_amqp_config(galaxy_config)
    with open(processing_config) as in_handle:
        config = yaml.load(in_handle)
    message_reader(analysis_handler(config), config["msg_tag"], amqp_config)

def copy_and_analyze(remote_info, config):
    """Remote copy an output directory, process it, and upload to Galaxy.
    """
    print remote_info

def message_reader(msg_handler, tag_name, config):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    conn = amqp.Connection(host=config['host'] + ":" + config['port'],
                           userid=config['userid'], password=config['password'],
                           virtual_host=config['virtual_host'], insist=False)
    chan = conn.channel()
    chan.queue_declare(queue=config['queue'], durable=True, exclusive=True,
            auto_delete=False)
    chan.exchange_declare(exchange=config['exchange'], type="direct", durable=True,
            auto_delete=False)
    chan.queue_bind(queue=config['queue'], exchange=config['exchange'],
                    routing_key=config['routing_key'])

    chan.basic_consume(queue=config['queue'], no_ack=True,
                       callback=msg_handler, consumer_tag=tag_name)
    while True:
        chan.wait()
    chan.basic_cancel(tag_name)
    chan.close()
    conn.close()

def analysis_handler(processing_config):
    tag_name = processing_config["msg_tag"]
    def receive_msg(msg):
        if msg.properties['application_headers'].get('msg_type') == tag_name:
            copy_and_analyze(json.loads(msg.body), processing_config)
    return receive_msg

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
