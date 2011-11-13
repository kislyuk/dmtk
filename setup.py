#!/usr/bin/env python

import os, sys, glob
from setuptools import setup, find_packages

setup(
    name='dmtk',
    version='0.1',
    description='DNA Modification Detection toolkit',
    author='Andrey Kislyuk',
    author_email='kislyuk@gmail.com',
    url='https://github.com/kislyuk/dmtk',
    license='GPL-3',
    package_dir = {'': 'lib'},
    packages = find_packages('lib'),
    scripts = glob.glob('scripts/*.py'),
    # fixme
    package_data = {'': ['R/*.R']},
    install_requires = ['numpy', 'scipy', 'h5py', 'matplotlib'],
)
