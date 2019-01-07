import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, debug
from dendropy import DnaCharacterMatrix
from copy import deepcopy


# run tiny_comb_... files first

workdir_ITS = "tests/data/PS_tiny_comb_its"
workdir_ETS = "tests/data/PS_tiny_comb_ets"
email = "mk@xy.zt"
percentage = 0.4

num_threads = 4  # number of threads to use, to make it run faster

pickle_fn = "final_ATT_checkpoint.p"

workdir_comb = "docs/example_scripts/output/nr_concat"
genelist = {"ITS": {"workdir": workdir_ITS, "pickle": pickle_fn}, 
            "ETS": {"workdir": workdir_ETS, "pickle": pickle_fn}
            }

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
                       email=email, num_threads=num_threads, percentage=percentage, user_concat_fn=None, backbone=None)