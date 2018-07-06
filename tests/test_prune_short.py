import os
import pickle
import sys

sys.stdout.write("\ntests prune_short\n")


try:
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

min_seqlen = 643

len_before = len(data_obj.tre.taxon_namespace)
data_obj.prune_short(min_seqlen)
len_after = len(data_obj.tre.taxon_namespace)

try:
	assert len_before > len_after
	sys.stdout.write("\nTEST passed: number of taxa in tre is shorter after pruning\n")
except:
    sys.stderr.write("\ntest failed\n")