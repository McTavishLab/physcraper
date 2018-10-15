import os
import sys
import pickle
#from physcraper import FilterBlast, ConfigObj, IdDicts
import physcraper
import physcraper.local_blast as local_blast


sys.stdout.write("\ntests read_local_blast\n")

# tests if I can read a local blast output file

workdir = "tests/output/test_read_local_blast"
configfi = "tests/data/test.config"
treshold = None
selectby = None
downtorank = None
absworkdir = os.path.abspath(workdir)

try:
    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape = physcraper.FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
filteredScrape.read_blast(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict()

# print("prepare test")
for taxonID in filteredScrape.sp_d:
    if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
        # print(taxonID)
        blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
        seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
        local_blast.write_blast_files(filteredScrape.workdir, taxonID, seq)
        blast_db = [item for item in filteredScrape.sp_seq_d[taxonID].keys()[1:] if type(item) == int]
        for blast_key in blast_db:
            seq = filteredScrape.sp_seq_d[taxonID][blast_key]
            local_blast.write_blast_files(filteredScrape.workdir, blast_key, seq, db=True, fn=str(taxonID))
        break

# test starts here:
blast_db = 1268580
blast_seq = 1268580
key = 1268580

local_blast.run_local_blast(filteredScrape.workdir, blast_seq, blast_db)
# print(filteredScrape.sp_seq_d.keys())
local_blast.read_local_blast(filteredScrape.workdir, filteredScrape.sp_seq_d[key], blast_db)

blast_out = "{}/blast/output_{}_tobeblasted.xml".format(workdir, key)

if os.path.exists(blast_out):
    with open(blast_out) as f:
        first_line = f.readline()
        try:
            assert len(first_line.strip()) != 0
            sys.stdout.write("\nTest passed!\n")
        # print("output file exists and is not empty, method can read the blast files and make an output file")
        except:
            sys.stderr.write("\ntest failed\n")
        # print(" output file of read_local_blast does not exist or is empty, method works not correctly")
