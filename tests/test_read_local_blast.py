import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can read a local blast output file


# I need to generate a FilterBlast object first
seqaln = "tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_read_local_blast"
configfi = "tests/data/aws.config"
id_to_spn = r"tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = None
selectby = None
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


if os.path.exists(otu_jsonfi):
    print("reload from file")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.isfile("{}/read_blast_test.p".format(workdir)): 
    filteredScrape = pickle.load(open("{}/read_blast_test.p".format(workdir),'rb'))
 
else:   
    conf = ConfigObj(configfi)
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                         mattype=mattype, 
                         workdir=workdir,
                         treefile=trfn,
                         schema_trf=schema_trf,
                         otu_json=otu_jsonfi,
                         ingroup_mrca=None)


    data_obj.prune_short()
    data_obj.dump()

    ids = IdDicts(conf, workdir=workdir)
    ids.dump()

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.run_blast()
    filteredScrape.read_blast()
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
    filteredScrape.dump('{}/read_blast_test.p'.format(workdir))

for taxonID in filteredScrape.sp_d:
	if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
	    print(taxonID)
	    blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
	    seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
	    filteredScrape.write_blast_files(taxonID, seq)

	    print(taxonID)
	    print(filteredScrape.sp_seq_d[taxonID].keys())
	    blast_db = [item for item in filteredScrape.sp_seq_d[taxonID].keys()[1:] if type(item) == int]
	    print(blast_db)
	    for blast_key in blast_db:
	    	seq = filteredScrape.sp_seq_d[taxonID][blast_key]

	    	filteredScrape.write_blast_files(blast_key, seq, db=True, fn=str(taxonID))
	    break





# test starts here:
blast_db = "Senecio_lagascanus"
blast_seq = "Senecio_lagascanus"
key = 'Senecio_lagascanus'

filteredScrape.run_local_blast(blast_seq, blast_db)
print(taxonID)

print(filteredScrape.sp_seq_d.keys())
filteredScrape.read_local_blast(filteredScrape.sp_seq_d[key], blast_db)


blast_out = "{}/blast/output_{}_tobeblasted.xml".format(workdir, key)

if os.path.exists(blast_out):
	with open(blast_out) as f:
		first_line = f.readline()
		try:
			assert len(first_line.strip()) != 0
			print("output file exists and is not empty, method can read the blast files and make an output file")
		except:
			print(" output file of read_local_blast does not exist or is empty, method works not correctly")