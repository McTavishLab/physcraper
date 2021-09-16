#!/usr/bin/env python

from setuptools import setup

setup(name='Physcraper'
      version='1.0.1',
      description='Physcraper',
      author='Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['physcraper'],
      scripts=['bin/physcraper_run.py',
               'bin/tree_comparison.py',
               'bin/find_trees.py',
               'bin/multi_loci.py'],
      long_description=(open('README.md').read()),
      url=("https://github.com/McTavishLab/physcraper"),
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
                        'nexson',
                        'numpy',
                        'wget']
     )
