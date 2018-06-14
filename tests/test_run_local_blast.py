import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can run a local blast query


# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_run_local_blast"
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = None
selectby = None
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

if os.path.isfile("{}/run_local_blast.p".format(workdir)): 
    filteredScrape = pickle.load(open("{}/run_local_blast.p".format(workdir),'rb'))
 
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
    filteredScrape.dump('{}/run_local_blast.p'.format(workdir))



# add var which you need to use for the test:
blast_db = "S_lagascanus"
blast_seq = "S_lagascanus"



filteredScrape.run_local_blast(blast_seq, blast_db)

blast_out = "{}/blast/output_S_lagascanus_tobeblasted.xml".format(workdir)
print(blast_out)
# xml_file = open(blast_out)
# print(xml_file)
if os.path.exists(blast_out):
    xml_file = open(blast_out)
    print(xml_file)
    print("succesfully read local blast output file")