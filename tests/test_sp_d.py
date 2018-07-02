import sys
import os
from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#

sys.stdout.write("\ntests sp_dict\n")

# tests if the building of sp_d and sp_seq_dict is working correclty
workdir = "tests/output/sp_d_test"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)
downtorank = None


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape =  FilterBlast(data_obj, ids)


try:
    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    filteredScrape.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    filteredScrape.read_blast(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]
    gi_data_otu_dict_added = []
    for v in filteredScrape.data.otu_dict.values():
        if '^ncbi:gi' in v:
            if (v['^physcraper:status'].split(' ')[0] not in filteredScrape.seq_filter):
                gi_data_otu_dict_added.append(v['^ncbi:gi'])
    gi_sp_d = []
    for key in filteredScrape.sp_d:
        v = filteredScrape.sp_d[key]  
        for v2 in v:
            if '^ncbi:gi' in v2:
                gi_sp_d.append(v2['^ncbi:gi'])
    user_data_otu_dict = []
    for v in filteredScrape.data.otu_dict.values():
        if '^user:TaxonName' in v:
            user_data_otu_dict.append(v['^user:TaxonName'])
    user_sp_d = []
    for v in filteredScrape.sp_d.values():
        for v2 in v:
            if '^user:TaxonName' in v2:
                user_sp_d.append(v2['^user:TaxonName'])
    assert sorted(gi_data_otu_dict_added) == sorted(gi_sp_d)
    assert sorted(user_data_otu_dict) == sorted(user_sp_d)
    # print("test passed: When generating the sp_d no entries which either have a gi or a userName are being lost.")
    sys.stdout.write("\ntest passed\n")
except:
    # print("failed: Test sp_d is not running correctly")
    sys.stderr.write("\ntest failed\n")

