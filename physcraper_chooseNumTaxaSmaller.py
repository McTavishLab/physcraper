from physcraper import wrappers_numTax
import os



#################################
seqaln =  "/home/martha/physcraper-git/physcraper/small_test_example/test.fas"
trfn= "/home/martha/physcraper-git/physcraper/small_test_example/test.tre"
id_to_spn = r"/home/martha/physcraper-git/physcraper/small_test_example/test_nicespl.csv"
workdir="numSpeciesSmall"
mattype="fasta"
schema_trf = "newick"
configfi = "example.config"
cwd = os.getcwd() 
treshhold=2


otu_json = wrappers_numTax.OtuJsonDict(id_to_spn, configfi)



wrappers_numTax.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
		treshhold,
                 otu_json,
                 configfi)
