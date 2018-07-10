import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


seqaln= "tests/data/tiny_test_example/test.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir="tests/output/owndata"
configfi = "example.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

"""Tests if your own input files will generate a data object of class AlignTreeTax
"""
sys.stdout.write("\nTesting 'generate_ATT_from_files'\n")
try:
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
								 ingroup_mrca=None)



	print(data_obj.gi_dict)
	print(data_obj.otu_dict)

	assert isinstance(data_obj, AlignTreeTax):
	sys.stdout.write("\nTest owndata.py passed\n")
	os.remove(otu_jsonfi)
#os.rmdir(workdir)
except:
	sys.stdout.write("\nTest owndata.py FAILED\n")
