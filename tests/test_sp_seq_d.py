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
workdir = "tests/output/sp_seq_d_test"
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

if os.path.isfile("{}/sp_seq_d_test.p".format(workdir)): 
    filteredScrape = pickle.load(open("{}/sp_seq_d_test.p".format(workdir),'rb'))
 
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
    filteredScrape.dump("{}/sp_seq_d_test.p".format(workdir))

filteredScrape.sp_dict(downtorank)

len_sp = 0
for k,v in filteredScrape.sp_d.iteritems():
    len_sp += len(v)

try:
    ## test begins
    gi_sp_d = []
    for key in filteredScrape.sp_d:
        v = filteredScrape.sp_d[key]  
        for v2 in v:
            if '^physcraper:status' in v2:
                not_added = ['deleted', 'subsequence,', 'not']
                if v2['^physcraper:status'].split(' ')[0] not in not_added: 
                    if '^ncbi:gi' in v2:
                        gi_sp_d.append(v2['^ncbi:gi'])

    user_sp_d = []
    count=0

    for v in filteredScrape.sp_d.values():
        for v2 in v:
            print(v2.keys())
            if '^physcraper:status' in v2 or u'^physcraper:status' in v2 :
                    not_added = ['deleted', 'subsequence,', 'not']
                    if v2['^physcraper:status'].split(' ')[0] not in not_added: 
                        if v2['^physcraper:last_blasted'] != '1800/01/01':
                            if '^user:TaxonName' in v2:
                                # print(v2['^ncbi:gi'])
                                user_sp_d.append(v2['^user:TaxonName'])
                                # print(v2['^user:TaxonName'])
                                count += 1
                            elif '^ot:ottTaxonName' in v2:
                                # print(v2['^ncbi:gi'])
                                user_sp_d.append(v2['^ot:ottTaxonName'])
                                # print(v2['^ot:ottTaxonName'])
                                count += 1


    filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
    print(filteredScrape.sp_seq_d.values())

    gi_sp_seq_d = []
    ott_sp_seq_d = []
    for v in filteredScrape.sp_seq_d.values():
        for k in v.keys():
            if type(k) == int:
                gi_sp_seq_d.append(k)
            if type(k) == str:
                # print(k)
                ott_sp_seq_d.append(k)

 
    assert len(ott_sp_seq_d) == len(user_sp_d)
    assert len(gi_sp_seq_d) == len(gi_sp_d)

    print("The length of the gi and user input names in sp_d and sp_seq_dict are the same, thus I expect the function to work correclty")

    print("Thus, transition from PhyscraperScrape to FilterBLAST should be working")

except:
    print("test sp_seq d not working!")
    # Problem might be related to different key encodings, but then I do not understand why it is not the same during the method of sp_seq_d