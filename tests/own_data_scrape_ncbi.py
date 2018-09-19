import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts

seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_ncbi_blast"
configfi = "tests/data/test.config"
# configfi = "tests/data/aws.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)


if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi, "w"))

wrappers.own_data_run(seqaln,
                      mattype,
                      trfn,
                      schema_trf,
                      workdir,
                      otu_jsonfi,
                      configfi)
