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
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None

print("trying to run sp_d")

try:
    if os.path.exists(otu_jsonfi):
        print("reload from file")
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
        if not os.path.exists(workdir):
           os.mkdir(workdir)
        json.dump(otu_json, open(otu_jsonfi,"w"))

    if os.path.isfile("{}/sp_d_test.p".format(workdir)): 
        filteredScrape = pickle.load(open("{}/sp_d_test.p".format(workdir),'rb'))
     
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
        filteredScrape.dump("{}/sp_d_test.p".format(workdir))

    filteredScrape.sp_dict(downtorank)
    print(filteredScrape.sp_d)
    # len_sp = 0
    # for k,v in filteredScrape.sp_d.iteritems():
    #     len_sp += len(v)

    # print("length sp_d:" , len_sp)


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

    print("When generating the sp_d no entries which either have a gi or a userName are being lost.")

except:
    print("Test sp_d is not running correctly")

