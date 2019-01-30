import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, IdDicts, ConfigObj

# tiny ets
seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"

workdir = "docs/example_scripts/output/tiny_comb_its"
configfi = "tests/data/localblast.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
threshold = 2
selectby = "blast"
downtorank = None

if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)

if os.path.exists(otu_jsonfi):
    print("load json")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = OtuJsonDict(id_to_spn, ids)
    json.dump(otu_json, open(otu_jsonfi, "w"))

wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         otu_jsonfi,
                         configfi,
                         selectby=selectby,
                         downtorank=downtorank)
