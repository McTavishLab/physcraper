import sys
import os
from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle
import shutil

# tests 
sys.stdout.write("\ntests loop_for_write_blast_files\n")


workdir = "tests/output/test_loop_for_write_blast_files"
configfi = "tests/data/test.config"
threshold = 2
selectby = "blast"
downtorank = None
absworkdir = os.path.abspath(workdir)

from pytest import mark

localblast = mark.localblast

@localblast
def test_loop_for_write_blast_files():
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.add_setting_to_self(downtorank, threshold)
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict()

    # this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
    # print("run loop which we want to test")
    for key in filteredScrape.sp_d:
        if len(filteredScrape.sp_d[key]) > threshold:
            count_dict = filteredScrape.count_num_seq(key)
            if key in filteredScrape.sp_seq_d.keys():
                seq_present = count_dict["seq_present"]
                query_count = count_dict["query_count"]
                if seq_present >= 1 and seq_present < threshold and count_dict["new_taxon"] is False and query_count != 0:
                    if query_count + seq_present > threshold:
                        taxonfn = filteredScrape.loop_for_write_blast_files(key)
                                    
# MAKE TEST FOR loop_for_write_blast_files

    print(filteredScrape.workdir)
    for key in filteredScrape.sp_d:
        count = 0
        count_int = 0
        count_gi_file = 0
        count_str_file = 0
        db = False
        blasted = False
        if len(filteredScrape.sp_d[key]) > threshold:
            for sp_keys in filteredScrape.sp_seq_d[key].keys():
                if isinstance(sp_keys, str):
                    count += 1
                if isinstance(sp_keys, unicode):
                    count += 1
                else:
                    count_int += 1
            folder = '{}/blast/'.format(filteredScrape.workdir)
            assert len(os.listdir(folder)) > 1, os.listdir(folder)
            for the_file in os.listdir(folder):
                # print(the_file)
                spn = int(the_file.split("_")[0])
                # print(spn)
                # spn = "_".join(the_file.split("_")[0])
                file_type = the_file.split("_")[1]
                assert spn in filteredScrape.sp_d.keys(), (spn, filteredScrape.sp_d.keys())

                print(spn, file_type)
                print(key)

                if spn == key:
                    if spn == key and file_type == "db": # 
                        db = True
                        f = open('{}/blast/{}'.format(filteredScrape.workdir, the_file))
                        for line in iter(f):
                            if line[0] == ">":
                                count_gi_file += 1
                    if spn == key and file_type == "tobeblasted":
                        blasted = True
                        count_str_file += 1
                    tested_both_files = 0
                    if blasted:
                        tested_both_files += 1
                        if count + count_int != threshold:
                            assert count_str_file == 1
                    if db:
                        tested_both_files += 1
                        if count + count_int != threshold:
                            assert count_gi_file == count_int
                
    assert tested_both_files >= 2

    shutil.rmtree('{}/blast/'.format(filteredScrape.workdir))

 