import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


#
seqaln= "small_test_example/test.fas"
mattype="fasta"
trfn= "small_test_example/test.tre"
schema_trf = "newick"
workdir="test_own_mini"
configfi = "tests/data/aws.config"
id_to_spn = r"small_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

cwd = os.getcwd()  

if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))



wrappers.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otu_jsonfi,
                 configfi)
