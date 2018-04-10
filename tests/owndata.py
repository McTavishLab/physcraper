import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


seqaln= "docs/owndata/senecio_its.fasta"
mattype="fasta"
trfn= "docs/owndata/its_new.tre"
schema_trf = "newick"
workdir="tests_owndata"
configfi = "example.config"
id_to_spn = r"docs/owndata/uniquetip_to_name_its.csv"
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
                             ingroup_mrca=None)


if isinstance(data_obj, AlignTreeTax):
    sys.stdout.write("{} is of type AlignTreeTax. Success, your input should be working.\n".format(data_obj))

os.remove(otu_jsonfi)
os.rmdir(workdir)
