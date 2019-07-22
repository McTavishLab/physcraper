import pickle
import sys
import os
from physcraper import ConfigObj, PhyscraperScrape, IdDicts

# Function we want to test is scrape.remove_identical_seqs()
# What are the inputs?
# physracper.scrape object, with new sequences read in.
# to make that we need: input data, idObject, and a configuration object.

# todo Make Sure
sys.stdout.write("Running test remove_identical_seqs\n\n")
workdir = "tests/data/tmp/owndata"
absworkdir = os.path.abspath(workdir)
conf = ConfigObj("tests/data/test.config", interactive=False)
conf.blast_loc='remote' #saves time over loading names and nodes, and they aren't used here


def test_remove_identical_seqs():
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir

    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    # print("start")
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.ids.otu_rank = {}
    scraper.config.gifilename = False
    scraper._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    #scraper.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    scraper.read_blast_wrapper(blast_dir=blast_dir)
    #print scraper.ncbi_mrca

    assert(len(scraper.new_seqs) == 0)
    assert(len(scraper.data.aln) == 5)
    assert len(scraper.new_seqs_otu_id) == 17
    #Now that we are pulling the full remote sequences, we don'thave any identical seuqnces in the test.

#TODO find an example where we do get identical sequences and need to discard them

    
    seqset = set()
    for otu in scraper.new_seqs_otu_id:
        seq = scraper.new_seqs_otu_id[otu]
        if seq in seqset:
            print otu
        seqset.add(seq)

#check that every new sequence is unique in the new seqs set, and is not a substring of another sequence.
    for otu in scraper.new_seqs_otu_id:
        qseq = scraper.new_seqs_otu_id[otu]
        count = 0
        for seq in seqset:
            if qseq in seq:
                count += 1
        assert count == 1


    for taxon in scraper.data.tre.taxon_namespace:
        assert(taxon.label in scraper.data.otu_dict)
        status = scraper.data.otu_dict[taxon.label].get(u'^physcraper:status')
        assert(status in ('original', 'query'))

#Now that we are pulling full sequences, length cutoff doensn't affect results.

#    # Second test checks that seq len prec is affecting results
#    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb')) #reload bc data object is mutable
#    data_obj.workdir = absworkdir
#    scraper2 = PhyscraperScrape(data_obj, ids)
#    scraper2.ids.otu_rank = {}

#    scraper2.config.gifilename = False
#    assert(len(scraper2.data.aln) == 5)
#    # scraper2.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
#    scraper2.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
#    scraper2.config.seq_len_perc = 2  # Change seq len percentage from default of 75%

#    assert(len(scraper2.new_seqs) == 40)
#    assert(len(scraper2.new_seqs_otu_id) == 0)

#    scraper2.remove_identical_seqs()
    # print(scraper2.data.otu_dict)
    # print(len(scraper.new_seqs_otu_id), 38)
    # print(len(scraper2.new_seqs_otu_id), 36)
#    assert(len(scraper2.new_seqs_otu_id) == 36)


  
test_remove_identical_seqs()