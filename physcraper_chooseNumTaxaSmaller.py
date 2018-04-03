from physcraper import wrappers_numTaxa
import os



#################################
seqaln =  "/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test.fas"
trfn= "/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test.tre"
id_to_spn = r"/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test_nicespl.csv"
workdir="numSpeciesSmaller"
mattype="fasta"
schema_trf = "newick"
configfi = "example.config"
cwd = os.getcwd() 
treshhold=2


otu_json = wrappers_numTaxa.OtuJsonDict(id_to_spn, configfi)



wrappers_numTaxa.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
		treshhold,
                 otu_json,
                 configfi)
