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
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


if os.path.exists(otu_jsonfi):
    print("reload from file")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.isfile("{}/select_seq_local_blast_test.p".format(workdir)): 
    print("reload to before test")
    filteredScrape = pickle.load(open("{}/select_seq_local_blast_test.p".format(workdir),'rb'))
 
else:   
    conf = ConfigObj(configfi)
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                         mattype=mattype, 
                         workdir=workdir,
                         treefile=trfn,
                         schema_trf=schema_trf,
                         otu_json=otu_jsonfi,
                         ingroup_mrca=None)

    data_obj.prune_short()
    data_obj.dump()

    ids = IdDicts(conf, workdir=workdir)
    ids.dump()

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.run_blast()
    filteredScrape.read_blast()
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
    filteredScrape.dump("{}/select_seq_local_blast_test.p".format(workdir))

##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
# filtered_seq = {}
print("run loop which we want to test")
for key in filteredScrape.sp_d:
    
    if len(filteredScrape.sp_d[key]) > treshold:
        print("filter number of sequences")
        print(key)
        count_dict = filteredScrape.count_num_seq(key)
        print(count_dict)
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
        for sp_keys in filteredScrape.sp_seq_d[key].keys():
            # print(type(sp_keys))
            if isinstance(sp_keys, str) == True:
                count += 1
            else:
                count_int +=1
        print(key)    
        folder = '{}/blast/'.format(filteredScrape.workdir)
        for the_file in os.listdir(folder):
            spn = "_".join(the_file.split("_")[1:-1])
            print("keys to compare")
            # print("_".join(key.split("_")[1:]))
            # print(spn)
            file_type = the_file.split("_")[-1]
            # print(file_type)
            # print("_".join(key.split("_")[1:]))
            if spn == "_".join(key.split("_")[1:]) and file_type == "db":
                print("db: names are equal")
                db = True
                f = open('{}/blast/{}'.format(filteredScrape.workdir, the_file))
                for line in iter(f):
                    # print(line)
                    if line[0] == ">":
                        count_gi_file += 1


                # go into the file and check for numbers of  >
            if spn == "_".join(key.split("_")[1:]) and file_type == "tobeblasted":
                print("tobeblasted: names are equal")
                blasted = True

                # print("something")
                count_str_file += 1
        if blasted:
            if count + count_int != treshold:
                print(count_str_file, count)

                assert count_str_file == count
        if db:
            
            if count + count_int != treshold:
                print(count_gi_file, count_int)
                assert count_gi_file == count_int


print("function loop_for_write_blast_files works")



