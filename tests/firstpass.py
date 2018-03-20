from physcraper.wrappers import standard_run
import pickle
import sys
import os


study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"
configfi = "lt.config"
workdir =  "november"

standard_run(study_id,
             tree_id,
             seqaln,
             mattype,
             workdir,
             configfi)
