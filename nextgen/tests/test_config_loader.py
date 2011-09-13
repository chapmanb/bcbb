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