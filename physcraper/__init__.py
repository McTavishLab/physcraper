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
from peyotl.api.phylesystem_api import PhylesystemAPI, APIWrapper
from peyotl.sugar import tree_of_life, taxomachine
from peyotl.nexson_syntax import (
    extract_tree,
    get_subtree_otus,
    extract_otu_nexson,
    PhyloSchema
)

import physcraper.AWSWWW as AWSWWW
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


def get_user_input():
    """Asks for yes or no user input.

    :return: user input
    """
    debug("get user input")
    is_valid = 0
    x = None
    while not is_valid:
        try:
            x = input("Please write either 'yes' or 'no': ")
            if x == "yes" or x == "no":
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError as e:
            print("'%s' is not a valid answer." % e.args[0].split(": ")[1])
    return x
