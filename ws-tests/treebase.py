from physcraper import ConfigObj, IdDicts,  PhyscraperScrape
from physcraper.opentree_helpers import get_dataset_from_treebase, generate_ATT_from_phylesystem, get_max_match_aln, get_tree_from_study

import pickle
import sys
import os

study_id = "ot_350"
tree_id = "Tr53297"
workdir = 'tests/output/treebase'
if not os.path.exists(workdir):
        os.makedirs(workdir)


alnfile = "{}/tb.aln".format(workdir)
aln_schema =  "nexus"

conf = ConfigObj()
tre, cite = get_tree_from_study(study_id, tree_id)
dataset = get_dataset_from_treebase(study_id)
aln = get_max_match_aln(tre, dataset)
aln.write(path=alnfile, schema = aln_schema)


data_obj = generate_ATT_from_phylesystem(alnfile=alnfile,
                                         aln_schema=aln_schema,
                                         workdir=workdir,
                                         configfile=conf,
                                         study_id=study_id,
                                         tree_id=tree_id)