from physcraper import wrappers
import os
import sys
import json


#################################

workdir="blast_numSpeciesTiny"
seqaln =  "tiny_test_example/test.fas"
trfn= "tiny_test_example/test.tre"
id_to_spn = r"tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

mattype="fasta"
schema_trf = "newick"
configfi = "tests/data/localblast.config"
cwd = os.getcwd() 
treshold=2
selectby="blast"
downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None



if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))


wrappers.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downtorank,
                 otu_jsonfi,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi)
