from physcraper import get_dataset_from_treebase, generate_ATT_from_phylesystem, ConfigObj, IdDicts,  PhyscraperScrape
import pickle
import sys
import os

study_id = "pg_873"
tree_id = "tree1679"
configfi = "tests/data/remotencbi.config"

conf = ConfigObj(configfi)

dataset = get_dataset_from_treebase(study_id,
                                phylesystem_loc='api')

aln = dataset.char_matrices[0]


data_obj = generate_ATT_from_phylesystem(aln=aln,
                                         workdir='tests/output/treebase',
                                         config_obj=conf,
                                         study_id=study_id,
                                         tree_id=tree_id)

