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


def test_sp_d():
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape =  FilterBlast(data_obj, ids)

    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
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
            v2 = filteredScrape.data.otu_dict[v2]
            if '^ncbi:gi' in v2:
                gi_sp_d.append(v2['^ncbi:gi'])
    user_data_otu_dict = []
    for v in filteredScrape.data.otu_dict.values():
        if '^user:TaxonName' in v:
            user_data_otu_dict.append(v['^user:TaxonName'])
    user_sp_d = []
    for v in filteredScrape.sp_d.values():
        for v2 in v:
            v2 = filteredScrape.data.otu_dict[v2]
            if '^user:TaxonName' in v2:
                user_sp_d.append(v2['^user:TaxonName'])
    assert sorted(gi_data_otu_dict_added) == sorted(gi_sp_d)
    assert sorted(user_data_otu_dict) == sorted(user_sp_d)
 