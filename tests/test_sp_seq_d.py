import sys
import os
from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle

sys.stdout.write("\ntests sp_seq_dict\n")


# tests if the building of sp_d and sp_seq_dict is working correclty
workdir = "tests/output/sp_seq_d_test"
configfi = "tests/data/test.config"
treshold = 2
selectby = "blast"
downtorank = None


absworkdir = os.path.abspath(workdir)
try:
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit(-1)
filteredScrape =  FilterBlast(data_obj, ids)


filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
# filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]

try:
    gi_sp_d = []
    for key in filteredScrape.sp_d:
        v = filteredScrape.sp_d[key]  
        for v2 in v:
            if '^physcraper:status' in v2:
                not_added = ['deleted', 'subsequence,', 'not']
                if v2['^physcraper:status'].split(' ')[0] not in not_added: 
                    if '^ncbi:gi' in v2:
                        gi_sp_d.append(v2['^ncbi:accession'])
    user_sp_d = []
    for v in filteredScrape.sp_d.values():
        for v2 in v:
            if '^physcraper:status' in v2 or u'^physcraper:status' in v2 :
                    if v2['^physcraper:status'].split(' ')[0] not in filteredScrape.seq_filter: 
                        if v2['^physcraper:last_blasted'] != '1800/01/01':
                            if '^user:TaxonName' in v2:
                                user_sp_d.append(v2['^user:TaxonName'])
                            elif '^ot:ottTaxonName' in v2:
                                user_sp_d.append(v2['^ot:ottTaxonName'])
    filteredScrape.make_sp_seq_dict()
    gi_sp_seq_d = []
    ott_sp_seq_d = []
    for v in filteredScrape.sp_seq_d.values():
        for k in v.keys():
            # print(k)
            if  len(k.split('.')) >=2:
            # if type(k) == int:
                gi_sp_seq_d.append(k)
            else:
            # if type(k) == str or type(k) == unicode:
                ott_sp_seq_d.append(k)
    # print(len(ott_sp_seq_d), len(user_sp_d), len(gi_sp_seq_d), len(gi_sp_d))
    assert len(ott_sp_seq_d) == len(user_sp_d)
    assert len(gi_sp_seq_d) == len(gi_sp_d)
    # print("The length of the gi and user input names in sp_d and sp_seq_dict are the same")
    sys.stdout.write("\ntest passed\n")
except:
    # print("FAILED:test sp_seq d not working!")
    sys.stderr.write("\ntest failed\n")
