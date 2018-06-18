import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can write proper labelled files


# I need to generate a FilterBlast object first
seqaln = "tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_write_label"
configfi = "tests/data/aws.config"
id_to_spn = r"tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = 'blast'
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None



otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
if not os.path.exists(workdir):
   os.mkdir(workdir)
json.dump(otu_json, open(otu_jsonfi,"w"))

  
conf = ConfigObj(configfi)
data_obj = generate_ATT_from_files(seqaln=seqaln, 
                     mattype=mattype, 
                     workdir=workdir,
                     treefile=trfn,
                     schema_trf=schema_trf,
                     otu_json=otu_jsonfi,
                     ingroup_mrca=None)

data_obj.write_labelled(label='user:TaxonName', treepath="labelled.tre", alnpath="labelled.fas")


if os.path.isfile('{}/labelled.tre'.format(workdir)) == True and os.path.isfile('{}/labelled.fas'.format(workdir)) == True:
	print('files written') 
else:
	print('did not write files')
	

