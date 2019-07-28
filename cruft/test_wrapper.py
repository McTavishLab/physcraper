import physcraper
import sys
import pickle
import os
from physcraper import ConfigObj, wrappers
from dendropy import DnaCharacterMatrix

#Use OpenTree phylesystem identifiers to get study and tree


def test_load_otol_data():
	study_id = "pg_873"
	tree_id = "tree1679"
	seqaln = "tests/data/minitest.fas"
	mattype = "fasta"
	workdir = "tests/output/opentree_unmappedtaxa"
	absworkdir = os.path.abspath(workdir)

	configfi = "tests/data/test.config"
	ingroup_mrca = None
	if not os.path.exists(workdir):
		os.mkdir(workdir)
	conf = ConfigObj(configfi)
	data_obj = wrappers.load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir)
	assert data_obj

def test_load_own_data():
	seqaln = "tests/data/tiny_test_example/test.fas"
	mattype = "fasta"
	trfn = "tests/data/tiny_test_example/test.tre"
	schema_trf = "newick"
	id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
	workdir = "tests/output/impls_mrcalist_local"
	configfi = "tests/data/test.config"
	ingroup_mrca = None

	if not os.path.exists(workdir):
		os.mkdir(workdir)

	conf = ConfigObj(configfi)

	ids = wrappers.load_ids_obj(conf, workdir)
	wrappers.make_otujsondict(id_to_spn, workdir, ids)

	data_obj = wrappers.load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
	assert data_obj



