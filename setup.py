#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='Physcraper',
      version='0.31',
      description='Physcraper',
      author='Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['physcraper'],
      scripts=['bin/physcraper_run.py',
               'bin/tree_comparison.py',
               'bin/find_trees.py'],
      include_package_data=True,
      install_requires=['argparse',
                        'biopython==1.76',
                        'configparser',
                        'DendroPy',
                        'DateTime',
                        'opentree',
                        'pandas',
                        'requests',
                        'sh',
                        'urllib3>=1.23',
                        'numpy',
                        'wget']
     )
