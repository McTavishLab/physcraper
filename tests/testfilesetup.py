import sys
import os
import json
import pickle
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, PhyscraperScrape
#


seqaln= "tests/data/tiny_test_example/test.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir="tests/output/owndata"
configfi = "tests/data/aws.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

"""Tests if your own input files will generate a data object of class AlignTreeTax
"""

if not os.path.exists("{}".format(workdir)):
	os.makedirs("{}".format(workdir))


otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
with open(otu_jsonfi,"w") as outfile:
    json.dump(otu_json, outfile)



data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             treefile=trfn,
                             schema_trf = schema_trf,
                             otu_json=otu_jsonfi,
                             email = "test@fake.edu",
                             ingroup_mrca=None)

data_obj.prune_short()
data_obj.dump(filename = "tests/data/tiny_dataobj.p")
conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)
scraper =  PhyscraperScrape(data_obj, ids, conf)
scraper.read_blast(blast_dir="tests/data/tiny_test_example/blast_files")
scraper.remove_identical_seqs()

pickle.dump(ids.gi_ncbi_dict, open("tests/data/tiny_gi_map.p", "wb" ))
