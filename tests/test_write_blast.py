#!/usr/bin/env python
import pickle
import sys
import os
from physcraper import ConfigObj, IdDicts
from physcraper import FilterBlast
import physcraper.local_blast as local_blast


sys.stdout.write("\ntests write_blast\n")

#set up test environment to test if method writes files for local blast 
workdir = "tests/output/test_write_local_blast_files"
configfi = "tests/data/test.config"

treshold = 2
selectby = "blast"
downtorank = None
absworkdir = os.path.abspath(workdir)

def test_write_blast():
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict()
    for taxonID in filteredScrape.sp_d:
        if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
            blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
            seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
            local_blast.write_filterblast_files(workdir, taxonID, seq)
            blast_db = filteredScrape.sp_seq_d[taxonID].keys()[1:]
            for blast_key in blast_db:
                seq = filteredScrape.sp_seq_d[taxonID][blast_key]
                local_blast.write_filterblast_files(workdir, blast_key, seq, db=True, fn=str(taxonID))
            break

    blast_file_blast = "{}/blast/{}_tobeblasted".format(workdir, taxonID)
    # print(blast_file_blast)
    blast_file_db = "{}/blast/{}_db".format(workdir, taxonID)
    # print(blast_file_db, blast_file_blast)
    if os.path.exists(blast_file_blast):
        with open(blast_file_blast) as f:
            first_line = f.readline()
            assert len(first_line.strip()) != 0
    if os.path.exists(blast_file_db):
        with open(blast_file_db) as f:
            first_line = f.readline()
            assert len(first_line.strip()) != 0
            
