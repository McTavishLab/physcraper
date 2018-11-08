import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, debug
from dendropy import DnaCharacterMatrix
from copy import deepcopy

def test():


    workdir_its = "./runs/tiny_comb_its"
    workdir_ets = "./runs/tiny_comb_ets"
    workdir_3 = "./runs/tiny_comb_ets"
    email = "martha.kandziora@yahoo.com"
    percentage = 0.4

    pickle_fn = "scrape_checkpoint.p"

    workdir_comb = "./tests/output/impl_concat"
    genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, 
                "ets": {"workdir": workdir_ets, "pickle": pickle_fn}
                }

    conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
                           email=email, percentage=percentage, user_concat_fn=None)
