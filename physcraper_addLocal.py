from physcraper import wrappers_numTax
import os



#################################
seqaln =  "/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test.fas"
trfn= "/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test.tre"
id_to_spn = r"/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test_nicespl.csv"
workdir="addLocal"
mattype="fasta"
schema_trf = "newick"
configfi = "example.config"
cwd = os.getcwd() 
treshold=10
selectby="blast" 
downto= None
add_local_seq = "/home/blubb/Documents/gitdata/physcraper/local_seqs"
id_to_spn_addseq = "/home/blubb/Documents/gitdata/physcraper/tipnTOspn_localAdd.csv"

otu_json = wrappers_numTax.OtuJsonDict(id_to_spn, configfi)

id_to_spn_addseq_json = wrappers_numTax.OtuJsonDict(id_to_spn_addseq, configfi)
# print(id_to_spn_addseq_json)

wrappers_numTax.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downto,
                 otu_json,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi)

#### is not adding sequences which are subsequences of sequences already present.