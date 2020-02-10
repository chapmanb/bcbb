#!/usr/bin/env python
"""Python setup file for Blue Collar Bioinformatics scripts and modules.
"""
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

__version__ = "Undefined"
for line in open('BCBio/GFF/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(name="bcbio-gff",
      version=__version__,
      author="Brad Chapman",
      author_email="chapmanb@50mail.com",
      license="Biopython License",
      description="Read and write Generic Feature Format (GFF) with Biopython integration.",
      url="https://github.com/chapmanb/bcbb/tree/master/gff",
      packages=find_packages(),
      install_requires=["six", "biopython"]
      )
