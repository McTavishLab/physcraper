import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can read a local blast output file


# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_read_local_blast"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = None
selectby = None
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)
filteredScrape._blasted = 1

blast_dir = "tests/data/precooked/fixed/tte_blast_files"
filteredScrape.read_blast(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)


print("prepare test")
# print(filteredScrape.sp_d)
for taxonID in filteredScrape.sp_d:
	if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
	    # print(taxonID)
	    blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
	    seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
	    filteredScrape.write_blast_files(taxonID, seq)
        # print("2nd taxonid")
        # print(taxonID)
        # print(filteredScrape.sp_seq_d[taxonID].keys())
        blast_db = [item for item in filteredScrape.sp_seq_d[taxonID].keys()[1:] if type(item) == int]
        # print(blast_db)
        for blast_key in blast_db:
	    	seq = filteredScrape.sp_seq_d[taxonID][blast_key]

	    	filteredScrape.write_blast_files(blast_key, seq, db=True, fn=str(taxonID))
        break




print("test begins")
# test starts here:
blast_db = "Senecio_lagascanus"
blast_seq = "Senecio_lagascanus"
key = 'Senecio_lagascanus'

filteredScrape.run_local_blast(blast_seq, blast_db)
# print(taxonID)

# print(filteredScrape.sp_seq_d.keys())
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