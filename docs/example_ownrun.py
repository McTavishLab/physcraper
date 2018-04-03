from physcraper import wrappers_nonstandard
import os

#
seqaln= "docs/owndata/senecio_its.fasta"
mattype="fasta"
trfn= "docs/owndata/its_new.tre"
schema_trf = "newick"
workdir="example_owndata_output_its"
configfi = "example.config"
id_to_spn = r"docs/owndata/uniquetip_to_name_its.csv"
cwd = os.getcwd()  



otu_json = wrappers_nonstandard.OtuJsonDict(id_to_spn, configfi)



wrappers_nonstandard.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otu_json,
                 configfi)
