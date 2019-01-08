import os
import sys
import pickle
#from physcraper import FilterBlast, ConfigObj, IdDicts
import physcraper
import physcraper.filter_by_local_blast as local_blast


sys.stdout.write("\ntests read_local_blast_query\n")

# tests if I can read a local blast output file

workdir = "tests/output/test_read_local_blast"
configfi = "tests/data/test.config"
treshold = None
selectby = None
downtorank = None
absworkdir = os.path.abspath(workdir)


def test_read_local_blast():
    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
    filteredScrape = physcraper.FilterBlast(data_obj, ids)
    filteredScrape._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    # filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    filteredScrape.sp_dict(downtorank)
    filteredScrape.make_sp_seq_dict()

    # print("prepare test")
    for taxonID in filteredScrape.sp_d:
        if len(filteredScrape.sp_seq_d[taxonID]) > treshold:
            # print(taxonID)
            blast_seq = filteredScrape.sp_seq_d[taxonID].keys()[0]
            seq = filteredScrape.sp_seq_d[taxonID][blast_seq]
            local_blast.write_filterblast_files(filteredScrape.workdir, taxonID, seq, fn=str(taxonID))
            # print(filteredScrape.sp_seq_d[taxonID].keys()[1:] )
            blast_db = [item for item in filteredScrape.sp_seq_d[taxonID].keys()[1:] if len(item.split(".")) >= 2]
            # print(blast_db)
            for blast_key in blast_db:
                seq = filteredScrape.sp_seq_d[taxonID][blast_key]
                local_blast.write_filterblast_files(filteredScrape.workdir, blast_key, seq, db=True, fn=str(taxonID))
            break

    print(taxonID) 
        
    blast_db = taxonID
    blast_seq = taxonID
    key = taxonID

    local_blast.run_filter_blast(filteredScrape.workdir, blast_seq, blast_db)
    local_blast.read_filter_blast(filteredScrape.workdir, filteredScrape.sp_seq_d[key], blast_db)

    blast_out = "{}/blast/output_{}_tobeblasted.xml".format(workdir, key)

    if os.path.exists(blast_out):
        with open(blast_out) as f:
            first_line = f.readline()
            assert len(first_line.strip()) != 0
          