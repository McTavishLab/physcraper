import sys
import os
import json
from physcraper import wrappers
import physcraper
from pytest import mark


def test_make_helper_classes():


	id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
	workdir = "tests/output/test_makeotu"
	configfi = "tests/data/localblast.config"


	if not os.path.exists(workdir):
		os.mkdir(workdir)
	conf = physcraper.ConfigObj(configfi)
	ids = wrappers.load_ids_obj(conf, workdir)
	assert os.path.exists("{}/id_pickle.p".format(workdir))

	wrappers.make_otujsondict(id_to_spn, workdir, ids)
	assert os.path.exists("{}/otu_dict.json".format(workdir))
	wrappers.make_otujsondict(id_to_spn, workdir, ids)

