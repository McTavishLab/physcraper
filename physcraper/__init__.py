#!/usr/bin/env python
"""Physcraper module"""

from __future__ import absolute_import

import sys
import re
import os
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import time
import csv
# from mpi4py import MPI
from past.builtins import xrange
from builtins import input
from copy import deepcopy
from ete2 import NCBITaxa


from Bio import Entrez

# extension functions
from . import concat  # is the local concat class
from . import ncbi_data_parser  # is the ncbi data parser class and associated functions
from . import filter_by_local_blast  # functions for the FilterBlast filtering
from . import opentree_helpers
from . import writeinfofiles


from physcraper.configobj import ConfigObj
from physcraper.ids import IdDicts
from physcraper.aligntreetax import AlignTreeTax
from physcraper.scrape import PhyscraperScrape
from physcraper.opentree_helpers import generate_ATT_from_phylesystem, OtuJsonDict, get_mrca_ott
from physcraper.aligntreetax import generate_ATT_from_files

if sys.version_info < (3,):
    from urllib2 import HTTPError
else:
    from urllib.error import HTTPError

_DEBUG = 1

_VERBOSE = 0

_DEBUG = 1
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)







# which python physcraper file do I use?
debug("Current --init-- version number: 12-17-2018.0")
debug(os.path.realpath(__file__))


