import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts

# define here your files
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
workdir = "tests/output/tiny_standard_own"
configfi = "tests/data/test.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

ingroup_mrca = None
shared_blast_folder = None

if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi, interactive=False)
ids = IdDicts(conf, workdir=workdir)


if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi, "w"))

# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.own_data_run(seqaln,
                      mattype,
                      trfn,
                      schema_trf,
                      workdir,
                      otu_jsonfi,
                      configfi,
		      ingroup_mrca=ingroup_mrca,
                      shared_blast_folder=shared_blast_folder)
