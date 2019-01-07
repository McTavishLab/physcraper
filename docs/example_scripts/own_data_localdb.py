from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts
import os
import json

#################################
seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
blacklist = None
workdir="tests/output/addLocal"

id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
otu_jsonfi_local = "{}/otu_dict_local.json".format(workdir)

configfi = "tests/data/localblast.config"
threshold=10
selectby="blast" 
downto= None
ingroup_mrca = None
add_unpubl_seq = "tests/data/local_seqs"
id_to_spn_addseq = "tests/data/tipnTOspn_localAdd.csv"


if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)


if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.exists(otu_jsonfi_local):
    print("load json local")
    otu_json_local = json.load(open(otu_jsonfi_local))
    print(otu_json_local)
else:
    otu_json_local = OtuJsonDict(id_to_spn_addseq, ids)
    json.dump(otu_json_local, open(otu_jsonfi_local,"w"))
    print(otu_json_local)

# print(id_to_spn_addseq_json)

wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         otu_jsonfi,
                         configfi,
                         selectby=selectby, 
                         downtorank=downto,
      			         ingroup_mrca=ingroup_mrca,
                         add_unpubl_seq=add_unpubl_seq,
                         id_to_spn_addseq_json=otu_json_local)
