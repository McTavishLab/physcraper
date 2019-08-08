import os
import json
import pickle
from physcraper import OtuJsonDict, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from physcraper import opentree_helpers
#


seqaln= "tests/data/tiny_test_example/test.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir="tests/data/tmp/owndata"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

"""Generates the files needed for the tests.
"""

if not os.path.exists("{}".format(workdir)):
	os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi, interactive=True)
ids = IdDicts(conf, workdir=workdir)

otu_json = OtuJsonDict(id_to_spn, ids)
with open(otu_jsonfi,"w") as outfile:
    json.dump(otu_json, outfile)


ottids = [otu_json[ite]['^ot:ottId'] for ite in otu_json]
mrca = opentree_helpers.get_mrca_ott(ottids)



data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             config_obj=conf,
                             treefile=trfn,
                             schema_trf = schema_trf,
                             otu_json=otu_jsonfi,
                             ingroup_mrca=mrca)

data_obj.prune_short()
data_obj.dump(filename = "tests/data/precooked/tiny_dataobj.p")

scraper =  PhyscraperScrape(data_obj, ids)

scraper._blasted = 1
scraper.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")

pickle.dump(ids.acc_ncbi_dict, open("tests/data/precooked/tiny_acc_map.p", "wb"))
# pickle.dump(scraper.acc_list_mrca, open("tests/data/precooked/acc_list_mrca.p", "wb"))

