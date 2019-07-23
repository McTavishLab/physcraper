import sys
import os
import pickle
import shutil
from physcraper import ConfigObj, IdDicts, PhyscraperScrape
from pytest import mark


from dendropy import Tree
import json

slow = mark.slow
localblast = mark.localblast


def new_test_generate_streamed_aln(Phy_obj):
    if Phy_obj.blacklist:
        Phy_obj.remove_blacklistitem()
    
    if len(Phy_obj.new_seqs) == 0 or len(Phy_obj.new_seqs_otu_id) == 0:
        Phy_obj.est_full_tree()
        Phy_obj.data.dump("{}/final_ATT_checkpoint.p".format(Phy_obj.workdir))
    elif len(Phy_obj.new_seqs) > 0:
        Phy_obj.data.write_files()  # should happen before aligning in case of pruning
        if len(Phy_obj.new_seqs_otu_id) > 0:  # TODO rename to something more intuitive
            Phy_obj.data.check_tre_in_aln()
            Phy_obj.write_query_seqs()
            Phy_obj.align_query_seqs()
            Phy_obj.place_query_seqs()
            Phy_obj.data.prune_short()
            Phy_obj.data.trim()
            Phy_obj.est_full_tree()
            Phy_obj.data.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(Phy_obj.workdir, Phy_obj.date),
                                     schema="newick",
                                     preserve_underscores=True,
                                     taxon_namespace=Phy_obj.data.aln.taxon_namespace)
            Phy_obj.data.write_files()
            
            Phy_obj.data.write_labelled(label='^physcraper:TaxonName', add_gb_id=True)
            Phy_obj.data.write_otus("otu_info", schema='table')
            Phy_obj.new_seqs = {}  # Wipe for next run
            Phy_obj.new_seqs_otu_id = {}
            Phy_obj.repeat = 1
        else:
            Phy_obj.est_full_tree()
            Phy_obj.data.dump("{}/final_ATT_checkpoint.p".format(Phy_obj.workdir))
    Phy_obj.reset_markers()
#    filter_by_local_blast.del_blastfiles(Phy_obj.workdir)  # delete local blast db
    Phy_obj.data.dump()
    json.dump(Phy_obj.data.otu_dict, open('{}/otu_dict.json'.format(Phy_obj.workdir), 'wb'))



@slow
@localblast
def test_blacklist():

    workdir = "tests/output/test_blacklist"
    configfi = "tests/data/test.config"

    # make one run without blacklist
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

    noblackScrape = PhyscraperScrape(data_obj, ids)
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
    new_test_generate_streamed_aln(noblackScrape)

    # one run with blacklist

    blacklist = ['JX895340.1']
    absworkdir = os.path.abspath(workdir)
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape = PhyscraperScrape(data_obj, ids)
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
    new_test_generate_streamed_aln(filteredScrape)


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
    new_test_generate_streamed_aln(noblackScrape)
    assert len_before - 1 == len(noblackScrape.data.tre.taxon_namespace)
 
