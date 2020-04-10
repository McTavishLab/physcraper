#!/usr/bin/env python
"""Physcraper module"""

from __future__ import absolute_import
import sys
import os

from physcraper.configobj import ConfigObj
from physcraper.ids import IdDicts
from physcraper.aligntreetax import AlignTreeTax
from physcraper.scrape import PhyscraperScrape
from physcraper.opentree_helpers import generate_ATT_from_phylesystem, OtuJsonDict, get_mrca_ott
from physcraper.aligntreetax import generate_ATT_from_files, generate_ATT_from_run
from physcraper.treetaxon import TreeTax

if sys.version_info < (3,):
    from urllib2 import HTTPError
else:
    from urllib.error import HTTPError

