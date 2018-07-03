import sys
import os
from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#

sys.stdout.write("\ntests select_seq_by_local_blast\n")

# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
workdir = "tests/output/test_select_seq_local_blast"
configfi = "tests/data/test.config"
treshold = 2
selectby = "blast"
downtorank = None

absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
filteredScrape.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
filteredScrape.read_blast(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)

##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
# print("start test")
count = 0
for giID in filteredScrape.sp_d:
    if len(filteredScrape.sp_d[giID]) > treshold:
        count_dict = filteredScrape.count_num_seq(giID)
        # print(count_dict)
        if count_dict["new_taxon"]:
            if count_dict["query_count"] < treshold:
                count += count_dict["query_count"]
            if count_dict["query_count"] > treshold:
                count += treshold
        if count_dict["new_taxon"] is False:
            if count_dict["seq_present"] < treshold:
                count += treshold-count_dict["seq_present"]
            if count_dict["seq_present"] > treshold:
                count += 0
        if giID in filteredScrape.sp_seq_d.keys():
            # print(giID in filteredScrape.sp_seq_d.keys())
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            # for item in filteredScrape.sp_d[giID]:
            if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
                if query_count + seq_present > treshold:
                    taxonfn = filteredScrape.loop_for_write_blast_files(giID, selectby)
                    for element in filteredScrape.sp_d[giID]:
                        if '^ot:ottTaxonName' in element:
                            blast_seq = "{}".format(element['^ot:ottTaxonName'])
                            blast_seq = blast_seq.replace(" ", "_")
                            blast_db = "{}".format(element['^ot:ottTaxonName'])
                            blast_db = blast_db.replace(" ", "_")
                    if filteredScrape.downtorank != None:
                        taxonfn = giID
                    filteredScrape.run_local_blast(taxonfn, taxonfn)
                    filteredScrape.select_seq_by_local_blast(filteredScrape.sp_seq_d[giID], taxonfn, treshold, seq_present)


try:
    assert count == len(filteredScrape.filtered_seq) and count>0
    sys.stdout.write("\ntest pass \n")
except:
    sys.stderr.write("\ntest failed\n")
