#!/usr/bin/env python
"""Setup script for distributed phylogenetic analysis by BLAST.
"""
from setuptools import setup, find_packages

setup(name = "bcbio-phyloblast",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "Distributed BLAST for gene phylogeny estimation",
      license = "BSD",
      url = "http://bcbio.wordpress.com",
      namespace_packages = ["bcbio"],
      packages = find_packages(),
      scripts = [],
      )
