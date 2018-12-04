import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, debug
from dendropy import DnaCharacterMatrix
from copy import deepcopy


# run tiny_comb_... files first

workdir_ITS = "./MS3_data/output/ITS_filter"
workdir_ETS = "./MS3_data/output/ETS_expand"
email = "mk@xy.zt"
percentage = 0.4

pickle_fn = "scrape_checkpoint.p"

workdir_comb = ".example/output/nr"
genelist = {"ITS": {"workdir": workdir_ITS, "pickle": pickle_fn}, 
            "ETS": {"workdir": workdir_ETS, "pickle": pickle_fn}
            }

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
                       email=email, percentage=percentage, user_concat_fn=None)

