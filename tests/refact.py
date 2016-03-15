from physcraper import StudyInfo, prune_short, config_obj, IdDicts
import pickle
import sys
import os


study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"
runname="fresh3"


configfi = "tests/local.config"

conf = config_obj(configfi)

info_obj = StudyInfo(seqaln,
                     mattype,
                     configfi,
                     study_id = study_id,
                     tree_id = tree_id)

data_obj = info_obj.generate_ATT()

prune_short(data_obj)
data_obj.write_files()
data_obj.write_labelled()


ids = IdDicts(conf.ott_ncbi, "tmp")