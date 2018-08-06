import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, debug
from dendropy import DnaCharacterMatrix
from copy import deepcopy#


#
workdir_its = "./tiny_comb_its"
workdir_ets = "./tiny_comb_ets"
workdir_3 = "./tiny_comb_ets"
email = "martha.kandziora@yahoo.com"
percentage = 0.4

pickle_fn = "scrape_checkpoint.p"


workdir_comb = "./tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, 
            "ets": {"workdir": workdir_ets, "pickle": pickle_fn}, 
            "3": {"workdir": workdir_3, "pickle": pickle_fn}}

##############
# print("{}/{}".format(workdir, pickle_fn))
# scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn),'rb'))
# print(scrape)


#############

# print(workdir_comb)

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb, 
						email=email, percentage = percentage, user_concat=None)

print(type(conc))


