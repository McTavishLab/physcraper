import os
import sys
import pickle
import filecmp

# tests if I can write proper labelled files

print("Running test_write_labelled\n\n")

expected_tree_path = "tests/data/expected_output/labelled.tre"
expected_aln_path = "tests/data/expected_output/labelled.fas"
expected_tree_path_ottid = "tests/data/expected_output/labelled_ottid.tre"
expected_aln_path_ottid = "tests/data/expected_output/labelled_ottid.fas"

try:
   data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
except:
   print("run 'python tests/testfilesetup.py' to setup data files for tests")

treepath = 'tests/data/tmp/labelled.tre'
alnpath = "tests/data/tmp/labelled.fas"

data_obj.write_labelled(label='^user:TaxonName', treepath=treepath, alnpath=alnpath, norepeats = False)

a =  os.path.isfile(treepath)
b = os.path.isfile(alnpath)
c = filecmp.cmp(treepath, expected_tree_path)
d = filecmp.cmp(alnpath, expected_aln_path)

treepath_ottid = 'tests/data/tmp/labelled_ottid.tre'
alnpath_ottid = "tests/data/tmp/labelled_ottid.fas"

data_obj.write_labelled(label='^ot:ottId', treepath=treepath_ottid, alnpath=alnpath_ottid, norepeats = False)

e = os.path.isfile(treepath_ottid)
f = os.path.isfile(alnpath_ottid)
g = filecmp.cmp(treepath_ottid, expected_tree_path_ottid)
h = filecmp.cmp(alnpath_ottid, expected_aln_path_ottid)

count = 0
if a*b*c*d*e*f*g*h:
	sys.stdout.write("Test write_lablled passed\n\n") 
else:
   sys.stderr.write("Test write_labelled FAILED\n\n")
   for var in [a,b,c,d,e,f,g,h]:
        if var != 1:
          	trials = ['a','b','c','d','e','f','g','h']
          	sys.stdout.write("{} is not true\n".format(trials[count]))
	count += 1
	

