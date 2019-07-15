import os
import sys
import pickle
import physcraper
import physcraper.filter_by_local_blast as local_blast
import physcraper.wrappers as wrappers
from physcraper.filterblast import FilterBlast

import pytest


def test_write_outputinfo():
    workdir = "tests/output/test_write_output_files"
    configfi = "tests/data/test.config"
    downtorank = None
    absworkdir = os.path.abspath(workdir)

    fn_otu = os.path.join(absworkdir, "otu_seq_info.csv")
    fn_sampling =  os.path.join(absworkdir,"taxon_sampling.csv")

    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    # filteredScrape.sp_dict(downtorank)
    # filteredScrape.make_sp_seq_dict()

    filteredScrape.generate_streamed_alignment()

    wrappers.write_out_files(filteredScrape, downtorank)

    with open(fn_otu) as fn:
        line = fn.readline()
        cnt = 1
        while cnt <= 5:
            line = fn.readline()
            cnt += 1
            assert type(line) == str       
            assert line.split(",") >= 2


    with open(fn_sampling) as fn:
        line = fn.readline()
        cnt = 1
        while cnt <= 5:
            line = fn.readline()
            cnt += 1
            assert type(line) == str
            assert line.split(",") >= 2
