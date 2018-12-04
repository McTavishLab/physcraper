import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts



seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "docs/example_scripts/output/own_standard_local"
configfi = "tests/data/localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)


ingroup_mrca = None
shared_blast_folder = None

if not os.path.exists("{}".format(workdir)):
	os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)

if os.path.exists(otu_jsonfi):
	print("load json")
else:
	otu_json = OtuJsonDict(id_to_spn, ids)
	json.dump(otu_json, open(otu_jsonfi, "w"))

# this function will keep all sequences found by blast which belong to the mrca, if you want to filter use filter_data_run()
wrappers.own_data_run(seqaln,
                  mattype,
                  trfn,
                  schema_trf,
                  workdir,
                  otu_jsonfi,
                  configfi,
                  ingroup_mrca=ingroup_mrca,
                  shared_blast_folder=shared_blast_folder)


