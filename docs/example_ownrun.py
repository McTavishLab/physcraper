from physcraper import wrappers
import os
import json

#
seqaln= "docs/owndata/senecio_its.fasta"
mattype="fasta"
trfn= "docs/owndata/its_new.tre"
schema_trf = "newick"
workdir="example_owndata_output_its"
configfi = "example.config"
id_to_spn = r"docs/owndata/uniquetip_to_name_its.csv"
otu_jsonfi = "example_owndata_output_its/otu_dict.json"

cwd = os.getcwd()  

if not os.path.exists(otu_jsonfi):
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    with open(otu_jsonfi,"w") as outfile:
        json.dump(otu_json, outfile)

wrappers.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otu_jsonfi,
                 configfi)
