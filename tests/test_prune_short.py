import os
import json
import pickle
import sys
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts 

## tests if prune_shorts, prunes seq from aln, tree and namespace
## not sure how to do that, as they need to be 


# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_read_local_blast"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = None
selectby = None
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()


#data_obj.prune_short()
# that part comes from prune_short method
min_seqlen = 643

print("test: prune short")

# print(len(data_obj.tre.taxon_namespace))
# print(len(data_obj.aln.taxon_namespace))
# print(data_obj.tre.as_string(schema='newick'))

len_before = len(data_obj.tre.taxon_namespace)
data_obj.prune_short(min_seqlen)

len_after = len(data_obj.tre.taxon_namespace)

# print(data_obj.tre.as_string(schema='newick'))
# print(len(data_obj.tre.taxon_namespace))
# print(len(data_obj.aln.taxon_namespace))

try:
	assert len_before != len_after
	print("number of taxa in tre is shorter after pruning")
except:
	print("test failed!")






