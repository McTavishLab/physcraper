from physcraper import wrappers_nonstandard
import os

#
seqaln= "Senecioneae/its_sl.fasta"
mattype="fasta"
trfn= "Senecioneae/its_sl.tre"
schema_trf = "newick"
workdir="Senecioneae_its_output"
configfi = "example.config"
id_to_spn = r"Senecioneae/nicespl.csv"
cwd = os.getcwd()  



otu_json = wrappers_nonstandard.OtuJsonDict(id_to_spn, configfi)



wrappers_nonstandard.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otu_json,
                 configfi)
