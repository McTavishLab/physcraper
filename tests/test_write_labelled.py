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

def test_write_labelled():
	data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))

	treepath = 'tests/data/tmp/labelled.tre'
	alnpath = "tests/data/tmp/labelled.fas"

	data_obj.write_labelled(label='^user:TaxonName', filename='labelled', direc='workdir', norepeats = False)

	a =  os.path.isfile(treepath)
	b = os.path.isfile(alnpath)
	c = filecmp.cmp(treepath, expected_tree_path)
	d = filecmp.cmp(alnpath, expected_aln_path)

	treepath_ottid = 'tests/data/tmp/labelled_ottid.tre'
	alnpath_ottid = "tests/data/tmp/labelled_ottid.fas"

	data_obj.write_labelled(label='^ot:ottId', filename='labelled_ottid', direc='tests/data/tmp/', norepeats = False)

	e = os.path.isfile(treepath_ottid)
	f = os.path.isfile(alnpath_ottid)
	g = filecmp.cmp(treepath_ottid, expected_tree_path_ottid)
	h = filecmp.cmp(alnpath_ottid, expected_aln_path_ottid)

	count = 0
	assert a*b*c*d*e*f*g*h == 1