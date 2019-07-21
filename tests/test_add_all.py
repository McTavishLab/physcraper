import sys
import os
import pickle  #
from physcraper import ConfigObj, IdDicts, PhyscraperScrape


workdir = "tests/output/add_all"
configfi = "tests/data/test.config"
threshold = 2
selectby = "blast"
downtorank = None
absworkdir = os.path.abspath(workdir)

def test_add_all():
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
   
    filteredScrape = PhyscraperScrape(data_obj, ids)
    filteredScrape._blasted = 1
    filteredScrape.threshold = threshold
    filteredScrape.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
    filteredScrape.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]
    filteredScrape.remove_identical_seqs()
    sp_d = filteredScrape.make_sp_dict(filteredScrape.new_seqs_otu_id)
    assert len(sp_d) == 7
    for taxon in sp_d:
        assert len(sp_d[taxon]) <= threshold