import sys
import os
import pickle  #
from physcraper import ConfigObj, IdDicts, FilterBlast

sys.stdout.write("\ntests add_all\n")

workdir = "tests/output/add_all"
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
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape = FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
filteredScrape.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict()
filteredScrape.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]
try:
    for key in filteredScrape.sp_d:
        if len(filteredScrape.sp_d[key]) <= treshold:
            filteredScrape.add_all(key)
    #############
    # print('test begins')
    treshold_undermin = 0
    for key in filteredScrape.sp_d:
        for key2 in filteredScrape.sp_d[key]:
            if len(filteredScrape.sp_d[key]) <= treshold:
                if '^physcraper:status' in key2:
                    if key2['^physcraper:status'].split(' ')[0] not in filteredScrape.seq_filter:
                        if key2['^physcraper:last_blasted'] == '1800/01/01':
                            treshold_undermin += 1
    add_all_thresholdmin = filteredScrape.filtered_seq
    assert treshold_undermin == len(add_all_thresholdmin)
    sys.stdout.write("\ntest passes\n")

except:
    sys.stderr.write("\ntest failed\n")
