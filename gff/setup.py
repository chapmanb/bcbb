#!/usr/bin/env python
"""Python setup file for Blue Collar Bioinformatics scripts and modules.
"""
from setuptools import setup, find_packages

__version__ = "Undefined"
for line in open('BCBio/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(name = "bcbio",
      version = __version__,
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "",
      packages = find_packages()
      #install_requires = ['SQLAlchemy>=0.4']
      )
