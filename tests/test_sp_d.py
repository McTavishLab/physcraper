import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast
import pickle#

# tests if the building of sp_d and sp_seq_dict is working correclty
#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/sp_d_test"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None

print("trying to run sp_d")

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

try:
    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    filteredScrape.read_blast(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()

    filteredScrape.sp_dict(downtorank)
    print(filteredScrape.sp_d)
    
    

    #attempt for test
    gi_data_otu_dict = []
    for v in filteredScrape.data.otu_dict.values():
        if '^ncbi:gi' in v:
            # print(v['^ncbi:gi'])
            gi_data_otu_dict.append(v['^ncbi:gi'])
    print(len(gi_data_otu_dict))


    gi_sp_d = []
    # print(filteredScrape.sp_d.keys())
    # print(filteredScrape.sp_d.values())
    for key in filteredScrape.sp_d:
        v = filteredScrape.sp_d[key]  

    # for v in filteredScrape.sp_d.values():
        # print(v)
        for v2 in v:
            if '^ncbi:gi' in v2:
                # print(v2['^ncbi:gi'])
                gi_sp_d.append(v2['^ncbi:gi'])
    print(len(gi_sp_d))


    user_data_otu_dict = []
    for v in filteredScrape.data.otu_dict.values():
        if '^user:TaxonName' in v:
            # print(v['^ncbi:gi'])
            user_data_otu_dict.append(v['^user:TaxonName'])
    print(len(user_data_otu_dict))

    user_sp_d = []
    # print(filteredScrape.sp_d.keys())
    # print(filteredScrape.sp_d.values())
    for v in filteredScrape.sp_d.values():
        # print(v)
        for v2 in v:
            if '^user:TaxonName' in v2:
                # print(v2['^ncbi:gi'])
                user_sp_d.append(v2['^user:TaxonName'])
    print(len(user_sp_d))



    assert sorted(gi_data_otu_dict) == sorted(gi_sp_d)
    assert sorted(user_data_otu_dict) == sorted(user_sp_d)
    # print(filteredScrape.sp_d)
    print("When generating the sp_d no entries which either have a gi or a userName are being lost.")

except:
    print("Test sp_d is not running correctly")

