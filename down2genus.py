from physcraper import wrappers_numTax
import os



#################################
seqaln= "/home/blubb/sync-TP-T470s/physcraper_testing/Senecioneae_input/its_sl.fasta"
mattype="fasta"
trfn= "/home/blubb/sync-TP-T470s/physcraper_testing/Senecioneae_input/its_sl.tre"
schema_trf = "newick"

id_to_spn = r"/home/blubb/sync-TP-T470s/physcraper_testing/Senecioneae_input/nicespl.csv"

workdir="blast2genus"
configfi = "example.config"
cwd = os.getcwd() 
treshhold=2
selectby="blast"
downtorank = "genus"
add_local_seq = None
id_to_spn_addseq_json = None




otu_json = wrappers_numTax.OtuJsonDict(id_to_spn, configfi)



wrappers_numTax.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
				treshhold,
				selectby,
				downtorank,
                 otu_json,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi)
