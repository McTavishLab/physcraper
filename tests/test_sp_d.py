import sys
import os
from physcraper import ConfigObj, IdDicts, PhyscraperScrape
import pickle#

sys.stdout.write("\ntests sp_dict\n")

# tests if the building of sp_d and sp_seq_dict is working correclty
workdir = "tests/output/sp_d_test"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)
downtorank = None


def test_sp_d():
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape =  PhyscraperScrape(data_obj, ids)

    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    filteredScrape.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]


