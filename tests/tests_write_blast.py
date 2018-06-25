#!/usr/bin/env python
import pickle
import sys
import os
import subprocess
import json
import csv
from ete2 import NCBITaxa
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape 
from physcraper import FilterBlast, debug #, Concat
from physcraper import wrappers
from dendropy import DnaCharacterMatrix


#set up test environment to test if method writes files for local blast 

# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_write_local_blast_files"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None

absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
filteredScrape.read_blast(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()


filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)


for taxonID in filteredScrape.sp_d:
	if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
	    blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
	    seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
	    filteredScrape.write_blast_files(taxonID, seq)

	    # print(taxonID)
	    # print(filteredScrape.sp_seq_d[taxonID].keys())
	    blast_db = filteredScrape.sp_seq_d[taxonID].keys()[1:]
	    print(blast_db)
	    for blast_key in blast_db:
	    	seq = filteredScrape.sp_seq_d[taxonID][blast_key]

	    	filteredScrape.write_blast_files(blast_key, seq, db=True, fn=str(taxonID))
	    break

print(taxonID)
blast_file_blast = "{}/blast/{}_tobeblasted".format(workdir, taxonID)
print(blast_file_blast)
blast_file_db = "{}/blast/{}_db".format(workdir, taxonID)
print(blast_file_db, blast_file_blast)
if os.path.exists(blast_file_blast):
	with open(blast_file_blast) as f:
		first_line = f.readline()
		try:
			assert len(first_line.strip()) != 0
			print("file exists and is not empty, method wrote the file for blast")
		except:
			print("file for blast does not exist or is empty, method works not correctly")


if os.path.exists(blast_file_db):
	with open(blast_file_db) as f:
		first_line = f.readline()
		try:
			assert len(first_line.strip()) != 0
			print("file exists and is not empty, method wrote the file for db")
		except:
			print("file for db does not exist or is empty, method works not correctly")