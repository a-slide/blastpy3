#!python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "blastpy3"
version = "0.3.0"
description = "Lightweight High level Python 3 API for NCBI BLAST"
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/a-slide/blastpy",
    author = 'Adrien Leger',
    author_email = 'aleg@ebi.ac.uk',
    license = 'GPLv3',
    python_requires ='>=3.5',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    packages = [name],
    package_dir = {name: name},
    install_requires = ["pyfaidx>=0.5.8"],
)
