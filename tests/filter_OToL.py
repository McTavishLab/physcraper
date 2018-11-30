import os
import json
from physcraper import wrappers, OtuJsonDict


# Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "tests/output/opentree"
configfi = "example.config"

threshold = 2
selectby = "blast"
downtorank = "species"

# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.filter_OTOL(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 configfi,
                 threshold,
                 selectby=selectby,
                 downtorank=downtorank)
