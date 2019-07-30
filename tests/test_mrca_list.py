# package import
import os
import json
import pickle
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts, PhyscraperScrape


from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 




def test_no_mrca():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = None
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=workdir)

    # print(ids.mrca_ott, ids.mrca_ncbi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = PhyscraperScrape(data_obj, ids, ingroup_mrca)
    filteredScrape.threshold = 5
    assert filteredScrape.mrca_ncbi == 18794

    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) == 24

