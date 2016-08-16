#!/usr/bin/env python

from setuptools import setup

setup(name='Physcraper',
      version='0.0',
      description='Physcraper',
      author='Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['physcraper'],
      dependency_links = ['https://github.com/OpenTreeOfLife/peyotl/tarball/physcraper#egg=peyotl-0.1.4devphyscraper'],
      install_requires=[
          'pickle',
          'dendropy',
          'configparser',
          'biopython',
          'urllib3',
          'peyotl==0.1.4devphyscraper'
      ]
     )
