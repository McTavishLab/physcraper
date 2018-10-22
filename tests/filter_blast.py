import sys
import os
import json
from physcraper import wrappers
#

seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "test_filter"
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
threshold = 2
selectby = "blast"

if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi, "w"))


wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         otu_jsonfi,
                         configfi,
                         selectby=selectby)
