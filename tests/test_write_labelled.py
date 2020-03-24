import os
import sys
import pickle
import filecmp

# tests if I can write proper labelled files

print("Running test_write_labelled\n\n")

expected_tree_path = "tests/data/expected_output/labelled.tre"
expected_aln_path = "tests/data/expected_output/labelled.fas"
"aesxpected_tree_path_ottid = "tests/data/expected_output/labelled_ottid.tre"
expected_aln_path_ottid = "tests/data/expected_output/labelled_ottid.fas"

def test_write_labelled():
	data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))

	treepath = 'tests/data/tmp/labelled.tre'
	alnpath = 'tests/data/tmp/labelled.fas'

	data_obj.write_labelled(label='^user:TaxonName', filename='labelled', direc='tests/data/tmp/', norepeats = False)

	assert os.path.isfile(treepath)
	assert os.path.isfile(alnpath)
	assert filecmp.cmp(treepath, expected_tree_path)
	assert filecmp.cmp(alnpath, expected_aln_path)
	os.remove(treepath)
	os.remove(alnpath)

	treepath_ottid = 'tests/data/tmp/labelled_ottid.tre'
	alnpath_ottid = 'tests/data/tmp/labelled_ottid.fas'

	data_obj.write_labelled(label='^ot:ottId', filename='labelled_ottid', direc='tests/data/tmp/', norepeats = False)

	assert os.path.isfile(treepath_ottid)
	assert os.path.isfile(alnpath_ottid)
	assert filecmp.cmp(treepath_ottid, expected_tree_path_ottid)
	assert filecmp.cmp(alnpath_ottid, expected_aln_path_ottid)
	os.remove(treepath_ottid)
	os.remove(alnpath_ottid)