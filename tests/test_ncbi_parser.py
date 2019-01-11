import pickle
import sys
import os
from physcraper import ConfigObj, PhyscraperScrape, IdDicts

import requests
import signal
from contextlib import contextmanager
from pytest import mark

slow = mark.slow


def test_ncbi_parser():

	workdir = "tests/output/test_run_raxml"
	absworkdir = os.path.abspath(workdir)
	conf = ConfigObj("tests/data/test.config", interactive=False)


	#load data 
	data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
	data_obj.workdir = absworkdir
	ids = IdDicts(conf, workdir=data_obj.workdir)
	ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

	taxid = ids.ncbi_parser.get_id_from_name("Crassocephalum vitellinum")
	print(taxid)
	
	rankidgenus = ids.ncbi_parser.get_downtorank_id(taxid, downtorank="genus")
	print(rankidgenus)

	mrcaid = ids.ncbi_parser.match_id_to_mrca(taxid, rankidgenus)
	print(mrcaid)

	rankid = ids.ncbi_parser.get_downtorank_id(taxid)
	print(rankid)

	

	rank = ids.ncbi_parser.get_rank(taxid)
	print(rank)

	synonym = ids.ncbi_parser.get_id_from_synonym("Elaps heterozonus")
	print(synonym)
	
	assert taxid == 1892268
	assert mrcaid == 189210
	assert rankid == 1892268
	assert rankidgenus == 189210

	assert rank == "species"
	assert synonym ==1117135


	