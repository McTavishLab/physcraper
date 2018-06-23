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
filteredScrape._blasted = 1
blast_dir = "tests/data/tiny_test_example/blast_files"
filteredScrape.read_blast(blast_dir="tests/data/precooked/fixed/tte_blast_files")
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)



try:
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
