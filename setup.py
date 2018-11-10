#!/usr/bin/env python

from distutils.core import setup

setup(name='Physcraper',
      version='0.0',
      description='Physcraper',
      author='Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['physcraper'],
      install_requires=[
          'dendropy==4.1.0',
          'configparser',
          'biopython',
          'urllib3',
          'peyotl',
      ]
     )
