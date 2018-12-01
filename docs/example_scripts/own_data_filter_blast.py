import sys
import os
import json
from physcraper import wrappers, OtuJsonDict
#


seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "docs/example_scripts/output/own_data_filter"
configfi = "tests/data/localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

threshold = 2
selectby = "blast"

ingroup_mrca = None  # must be OToL ID
shared_blast_folder = None # location to share blast runs across runs, see documentation
downtorank = None
blacklist = None
add_unpubl_seq = None
id_to_spn_addseq_json = None

if not os.path.exists(workdir):
    os.mkdir(workdir)


if os.path.exists(otu_jsonfi):
	otu_json = json.load(open(otu_jsonfi))
else:
	otu_json = OtuJsonDict(id_to_spn, configfi)
	json.dump(otu_json, open(otu_jsonfi, "w"))



wrappers.filter_data_run(seqaln,
                     mattype,
                     trfn,
                     schema_trf,
                     workdir,
                     threshold,
                     otu_jsonfi,
                     configfi,
                     downtorank=downtorank,
                     selectby=selectby,
					 blacklist=blacklist,
                     add_unpubl_seq=add_unpubl_seq,
                     id_to_spn_addseq_json=id_to_spn_addseq_json,
                     ingroup_mrca=ingroup_mrca,
                     shared_blast_folder=shared_blast_folder
                     )

