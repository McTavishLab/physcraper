from physcraper import generate_ATTs_from_treebase, generate_ATT_from_files, prune_short, ConfigObj, IdDicts,  PhyscraperScrape
import pickle
import sys
import os


study_id = "pg_873"
tree_id = "tree1679"

configfi = "tests/local.config"

conf = ConfigObj(configfi)

dataset = get_dataset_from_treebase(study_id,
                                phylesystem_loc='api'):


aln = dataset[0]

ATTs = generate_ATTs_from_treebase(aln=aln,
                                  workdir='treebase',
                                  study_id=study_id,
                                  tree_id=tree_id,
                                  phylesystem_loc='local'):