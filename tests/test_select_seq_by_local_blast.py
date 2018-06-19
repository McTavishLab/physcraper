import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast
import pickle#

# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_select_seq_local_blast"
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
print("start test")
count = 0
for giID in filteredScrape.sp_d:
    if len(filteredScrape.sp_d[giID]) > treshold:
        count_dict = filteredScrape.count_num_seq(giID)
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
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            for item in filteredScrape.sp_d[giID]:
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
                    


print(count)
print(len(filteredScrape.filtered_seq))

try:
    assert count == len(filteredScrape.filtered_seq)
    print("test passed")
except:
    print("test failed")