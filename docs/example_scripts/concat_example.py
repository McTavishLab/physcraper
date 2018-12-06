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

pickle_fn = "scrape_checkpoint.p"

workdir_comb = "docs/example_scripts/output/nr_concat"
genelist = {"ITS": {"workdir": workdir_ITS, "pickle": pickle_fn}, 
            "ETS": {"workdir": workdir_ETS, "pickle": pickle_fn}
            }

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
                       email=email, percentage=percentage, user_concat_fn=None)

