#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='Physcraper',
      version='0.1',
      description='Physcraper',
      author='Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['physcraper'],
      scripts=['bin/opentree_scrape.py'],
      install_requires=[]
     )
