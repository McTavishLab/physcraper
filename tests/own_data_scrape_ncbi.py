import sys
import os
import json
import shutil
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


#
seqaln= "tests/data/tiny_test_example/test.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir="tests/test_ncbi_blast"
configfi = "tests/data/ncbiblast.config"
# configfi = "tests/data/aws.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold=2
selectby="blast"
downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None

cwd = os.getcwd()  

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

shutil.rmtree(workdir)