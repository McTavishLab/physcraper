import os
import json
from physcraper import wrappers, OtuJsonDict

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"
workdir="tests/output/opentree"
configfi = "example.config"


otu_jsonfi = "{}/otu_dict.json".format(workdir)
threshold = 2
selectby = "blast"
downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None

print(workdir)
if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))


# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.filter_OTOL(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi,
                 treshold,
                 selectby,
                 downtorank,
                 spInfoDict,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 blacklist)
