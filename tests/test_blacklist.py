import sys
import os
import pickle
import shutil
import pytest
from physcraper import ConfigObj, IdDicts, FilterBlast, debug

from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
slow = mark.slow


@slow
def test_blacklist():

    workdir = "tests/output/test_blacklist"
    configfi = "tests/data/test.config"

    # make one run without blacklist
    debug("run without blacklist")
    blacklist = None
    noblack = os.path.join(workdir, "noblacklist")
    absworkdir = os.path.abspath(noblack)
    if not os.path.exists(os.path.join(absworkdir, "current_blast_run/")):
        os.makedirs(os.path.join(absworkdir, "current_blast_run/"))


    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    noblackScrape = FilterBlast(data_obj, ids)
    noblackScrape._blasted = 1
    src = "tests/data/precooked/fixed/tte_blast_files"
    src_files = os.listdir(src)
    for file_name in src_files:
        dest = os.path.join(absworkdir, "current_blast_run/")
        # print(dest)
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    noblackScrape.read_blast_wrapper()
    noblackScrape.remove_identical_seqs()
    noblackScrape.generate_streamed_alignment()

    # one run with blacklist
    debug("run with blacklist")

    blacklist = ['JX895340.1']
    absworkdir = os.path.abspath(workdir)
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.blacklist = blacklist
    filteredScrape._blasted = 1
    if not os.path.exists(os.path.join(absworkdir, "current_blast_run/")):
        os.makedirs(os.path.join(absworkdir, "current_blast_run/"))
    src = "tests/data/precooked/fixed/tte_blast_files"
    src_files = os.listdir(src)
    for file_name in src_files:
        dest = os.path.join(absworkdir, "current_blast_run/")
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper()
    filteredScrape.remove_identical_seqs()
    filteredScrape.generate_streamed_alignment()


    print("RUN TESTS!")
    gi_l = []
    gi_l_2 = []
    for tax in filteredScrape.data.tre.taxon_namespace:
        gi_id = filteredScrape.data.otu_dict[tax.label].get("^ncbi:accession")
        gi_l.append(gi_id)
    print(gi_l)
    for tax in noblackScrape.data.tre.taxon_namespace:
        # print(filteredScrape.data.otu_dict[tax.label])
        gi_id = noblackScrape.data.otu_dict[tax.label].get("^ncbi:accession")
        gi_l_2.append(gi_id)
    print(gi_l_2)
    for item in blacklist:
        assert item not in gi_l
        print("RUN TESTS2!")
        assert item in gi_l_2
    
        #     # print("seq was not added in blacklist run")
        #     print("inbetween step works")
# test if it removes blacklist gi from already added aln:
    print("run with later blacklist")

    # else:
    #     print("blacklist gi was added in previous run")
    # print("now we want to remove it.")
    len_before = (len(noblackScrape.data.tre.taxon_namespace))
    noblackScrape.blacklist = blacklist
    noblackScrape.generate_streamed_alignment()
    assert len_before - 1 == len(noblackScrape.data.tre.taxon_namespace)
 