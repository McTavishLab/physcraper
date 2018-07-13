import os
import json
from physcraper import wrappers, OtuJsonDict
#


#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
workdir="tiny_run_july"
configfi = "tests/data/blubb_localblast_highhits.config"


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
wrappers.own_data_run(seqaln,
                     mattype,
                     trfn,
                     schema_trf,
                     workdir,
                     downtorank,
                     otu_jsonfi,
                     configfi)

