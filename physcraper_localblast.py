from physcraper import wrappers_numTax
import os



#################################
seqaln =  "/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test.fas"
trfn= "/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test.tre"
id_to_spn = r"/home/blubb/Documents/gitdata/physcraper/smaller_test_example/test_nicespl.csv"
workdir="localblast_numSpeciesSmaller"
mattype="fasta"
schema_trf = "newick"
configfi = "example.config"
cwd = os.getcwd() 
treshhold=2
selectby="blast"


otu_json = wrappers_numTax.OtuJsonDict(id_to_spn, configfi)



wrappers_numTax.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
				treshhold,
				selectby,
                 otu_json,
                 configfi)
