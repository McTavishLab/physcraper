from physcraper import wrappers_nonstandard
import os

#
seqaln= "/home/blubb/sync-TP-T470s/physcraper_testing/concatenate_aln/labelled_ets.fas"
mattype="fasta"
trfn= "/home/blubb/sync-TP-T470s/physcraper_testing/concatenate_aln/labelled_ets.tre"
schema_trf = "newick"
workdir= "/home/blubb/sync-TP-T470s/physcraper_testing/senecio_out_ets"
configfi = "example.config"
id_to_spn = r"/home/blubb/sync-TP-T470s/physcraper_testing/concatenate_aln/nicespl.csv"
cwd = os.getcwd()  



otu_json = wrappers_nonstandard.OtuJsonDict(id_to_spn, configfi)



wrappers_nonstandard.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otu_json,
                 configfi)
