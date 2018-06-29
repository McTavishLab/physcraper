import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast
import pickle#

# tests loop_for_write_blast_files
#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_loop_for_write_blast_files"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)
filteredScrape.read_blast(blast_dir="tests/data/precooked/fixed/tte_blast_files")
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)

##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
# filtered_seq = {}
print("run loop which we want to test")
for key in filteredScrape.sp_d:
    
    if len(filteredScrape.sp_d[key]) > treshold:
        # print("filter number of sequences")
        # print(key)
        count_dict = filteredScrape.count_num_seq(key)
        # print(count_dict)
        if key in filteredScrape.sp_seq_d.keys():
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
                    if query_count + seq_present > treshold:
                        taxonfn = filteredScrape.loop_for_write_blast_files(key, selectby)
                                


####MAKE TEST FOR loop_for_write_blast_files
print("run the test")
for key in filteredScrape.sp_d:
    count = 0
    count_int = 0
    count_gi_file = 0
    count_str_file = 0
    db = False
    blasted = False
    if len(filteredScrape.sp_d[key]) > treshold:
        # print(filteredScrape.sp_seq_d[key])
        for sp_keys in filteredScrape.sp_seq_d[key].keys():
            # print(type(sp_keys))
            if isinstance(sp_keys, str) == True:
                count += 1
                # blasted = True
            if isinstance(sp_keys, unicode) == True:
                count += 1
                # blasted = True
            else:
                count_int +=1
                # db = True

        # print(key)    
        folder = '{}/blast/'.format(filteredScrape.workdir)
        for the_file in os.listdir(folder):
            # print(the_file)
            spn = "_".join(the_file.split("_")[1:-1])
            # print("keys to compare")
            # print("_".join(key.split("_")[1:]))
            # print(spn)
            file_type = the_file.split("_")[-1]
            # print(file_type)
            # print("_".join(key.split("_")[1:]))
            if spn == "_".join(key.split("_")[1:]) and file_type == "db":
                # print("db: names are equal")
                db = True
                f = open('{}/blast/{}'.format(filteredScrape.workdir, the_file))
                for line in iter(f):
                    # print(line)
                    if line[0] == ">":
                        count_gi_file += 1


                # go into the file and check for numbers of  >
            if spn == "_".join(key.split("_")[1:]) and file_type == "tobeblasted":
                # print("tobeblasted: names are equal")
                blasted = True

                # print("something")
                count_str_file += 1
        # print(count, count_int, count_str_file, treshold)
        # print(blasted)
        

        if blasted:
            if count + count_int != treshold:
                print(count_str_file, count)

                assert count_str_file == count
        if db:
            # print(count, count_int, treshold)
            if count + count_int != treshold:
                print(count_gi_file, count_int)
                assert count_gi_file == count_int


print("function loop_for_write_blast_files works")



