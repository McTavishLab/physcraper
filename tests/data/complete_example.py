# package import
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts

# define here your files
seqaln = "tests/data/Senecioneae_input/its_sl.fasta"
mattype = "fasta"
trfn = "tests/data/Senecioneae_input/its_sl.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/Senecioneae_input/nicespl.csv"
workdir="runs_for_MS/Senecioneae_filter4"
configfi = "tests/data/test.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
shared_blast_folder = "/home/mkandziora/physcraper/shared_runs/"


# change to your filtering criteria
threshold = 4
selectby = "blast"
downtorank = "genus"
add_unpubl_seq = None  # path to file
id_to_spn_addseq = None # path to file
otu_jsonfi_local = "{}/otu_dict_unpubl.json".format(workdir)

blacklist = None  # is a list with gi numbers
ingroup_mrca = 557768

# setup the run
if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)


if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi,"w"))

if add_unpubl_seq is not None:
    if os.path.exists(otu_jsonfi_local):
        print("load json local")
        otu_json_local = json.load(open(otu_jsonfi_local))
        print(otu_json_local)
    else:
        otu_json_local = OtuJsonDict(id_to_spn_addseq, ids)
        json.dump(otu_json_local, open(otu_jsonfi_local,"w"))
        print(otu_json_local)


# select a wrapper function, depending on what you want to do, see short tutorial:
wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         otu_jsonfi,
                         configfi,
                         selectby=selectby,  # default is "blast", "length" is the other option
                         downtorank=downtorank,  # default is "species", can be any hierachy which exists with name in the taxonomy browser of ncbi
                         blacklist=blacklist,
                         add_unpubl_seq=add_unpubl_seq,
                         id_to_spn_addseq_json=id_to_spn_addseq,
                         ingroup_mrca=ingroup_mrca,
                         shared_blast_folder=shared_blast_folder)




