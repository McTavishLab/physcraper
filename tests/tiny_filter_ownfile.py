# package import
import os
import json
from physcraper import wrappers, OtuJsonDict

# define here your files
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
workdir="tiny_run_july"
configfi = "tests/data/blubb_localblast.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

# change to your filtering criteria
treshold = 2
selectby = "blast"
downtorank = "species"
add_local_seq = None  # under development
id_to_spn_addseq_json = None  # under development
blacklist = None
# setup the run
if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))


# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         treshold,
                         selectby,
                         downtorank,
                         otu_jsonfi,
                         blacklist,
			 add_local_seq,
                         id_to_spn_addseq_json,
                         configfi)

