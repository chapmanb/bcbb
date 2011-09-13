"""Tests the bcbio.pipeline.config_loader module.
"""

import os
import sys

sys.path.append(os.path.realpath(".."))
from bcbio.pipeline.config_loader import load_config

def test_loading():
	"""Test loading a given file.
	"""
	config = load_config("data/loading_test/variables.yaml")
	assert(isinstance(config, dict))

def test_variable_expansion():
	"""Test expanding the environment variables in the
	test yaml.
	"""
	config = load_config("data/loading_test/variables.yaml")

	try:
		for variable, value in config.items():
			assert (os.environ[variable] == value
			), "The strings %s and %s doesn't match (variable %s)" % (
			os.environ[variable], value, variable)
	
	# When the key isnn't in os.environ
	except KeyError as e:
		for variable, value in config[e.message].items():
			assert (os.environ[variable] == value
			), "The strings %s and %s doesn't match (variable %s)" % (
			os.environ[variable], value, variable)
