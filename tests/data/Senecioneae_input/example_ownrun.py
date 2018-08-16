from physcraper import wrappers, ConfigObj, IdDicts, OtuJsonDict
import os

#
seqaln= "tests/data/Senecioneae_input/its_sl.fasta"
mattype="fasta"
trfn= "tests/data/Senecioneae_input/its_sl.tre"
schema_trf = "newick"
workdir="test_iddict"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/Senecioneae_input/nicespl.csv"
cwd = os.getcwd()  
otu_jsonfi = "{}/otu_dict.json".format(workdir)


treshold = 2
selectby = "blast"
downtorank = None
blacklist = None
# downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None



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
