import sys
import os
import json
import pickle#
from physcraper import debug, wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast


print('tests add_all')
# tests if the building of select_seq_from local blast is selecting the right amount of species
#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/add_all"
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None

# if os.path.exists(otu_jsonfi):
#     print("reload from file")
#     otu_json = json.load(open(otu_jsonfi))
# else:
#     otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
#     if not os.path.exists(workdir):
#        os.mkdir(workdir)
#     json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.isfile("{}/test_add_all.p".format(workdir)): 
    print("reload to before test")
    filteredScrape = pickle.load(open("{}/test_add_all.p".format(workdir),'rb'))
 
else:   
    data_obj = pickle.load(open("tests/data/tiny_dataobj.p", 'rb'))
    data_obj.workdir = workdir
    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/tiny_gi_map.p", "rb" ))
    filteredScrape =  FilterBlast(data_obj, ids)
    filteredScrape._blasted = 1
    blast_dir = "tests/data/tiny_test_example/blast_files"
    filteredScrape.read_blast(blast_dir=blast_dir)

    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
    filteredScrape.dump("{}/test_add_all.p".format(workdir))




try:
    # data_obj = pickle.load(open("tests/data/tiny_dataobj.p", 'rb'))
    # conf = ConfigObj(configfi)
    # ids = IdDicts(conf, workdir=data_obj.workdir)
    # ids.gi_ncbi_dict = pickle.load(open("tests/data/tiny_gi_map.p", "rb" ))
    # filteredScrape =  FilterBlast(data_obj, ids)
    # filteredScrape._blasted = 1
    # blast_dir = "tests/data/tiny_test_example/blast_files"
    # filteredScrape.read_blast(blast_dir=blast_dir)

    # filteredScrape.sp_dict(downtorank)
    # filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
    # filteredScrape.dump("test_add_all.p")


    for key in filteredScrape.sp_d:
        if len(filteredScrape.sp_d[key]) <= treshold:
            filteredScrape.add_all(key)

    #############
    ### make actual test
    print('test begins')
    treshold_undermin = 0
    for key in filteredScrape.sp_d:
        for key2 in filteredScrape.sp_d[key]:
            if len(filteredScrape.sp_d[key]) <= treshold:
                if '^physcraper:status' in key2:
                    not_to_add = ['deleted', 'subsequence,', 'not']

                    if key2['^physcraper:status'].split(' ')[0] not in not_to_add: 
                        if key2['^physcraper:last_blasted'] == '1800/01/01':
                            treshold_undermin += 1
    add_all_thresholdmin = filteredScrape.filtered_seq
    try:
        assert treshold_undermin == len(add_all_thresholdmin) 
        print("test for add_all passes")
    except:
        print("test for add_all does not pass")

except:
    print("try failed...")
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))


    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    with open(otu_jsonfi,"w") as outfile:
        json.dump(otu_json, outfile)


    conf = ConfigObj(configfi)
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                 mattype=mattype, 
                                 workdir=workdir,
                                 treefile=trfn,
                                 schema_trf = schema_trf,
                                 otu_json=otu_jsonfi,
#                                 email = conf.email,
                                 ingroup_mrca=None)

    data_obj.prune_short()
    data_obj.dump(filename = "tests/data/tiny_dataobj.p")

    ids = IdDicts(conf, workdir=workdir)
    ids.dump()
    filteredScrape =  FilterBlast(data_obj, ids)
    filteredScrape.read_blast(blast_dir="tests/data/tiny_test_example/blast_files")
    filteredScrape.remove_identical_seqs()

    pickle.dump(ids.gi_ncbi_dict, open("tests/data/tiny_gi_map.p", "wb" ))
    print('rerun the test, the files were not present')
    