import pickle
import sys
import os
from physcraper import ConfigObj, PhyscraperScrape, IdDicts

#Function we want to test is scrape.remove_identical_seqs()
#What are the inputs?  
#physracper.scrape object, with new sequences read in.
#to make that we need: input data, idObject, and a configuration object.

##todo Make Sure 

workdir="tests/output/owndata"
absworkdir = os.path.abspath(workdir)

sys.stdout.write("\nTesting 'remove_identical_seqs'\n")
try:
    data_obj = pickle.load(open("tests/data/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    conf = ConfigObj("tests/data/aws.config")
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/tiny_gi_map.p", "rb" ))
    scraper =  PhyscraperScrape(data_obj, ids)
    scraper._blasted = 1
    blast_dir = "tests/data/tiny_test_example/blast_files"
    scraper.read_blast(blast_dir=blast_dir)

    assert len(scraper.new_seqs) == 40
    assert len(scraper.data.aln) == 5
    assert len(scraper.new_seqs_otu_id) == 0
    
    scraper.remove_identical_seqs()

    assert len(scraper.new_seqs) == 40
    assert len(scraper.data.aln) == 5
    assert len(scraper.new_seqs_otu_id) == 37


#Second test checks that seq len prec is affecting results
    data_obj = pickle.load(open("tests/data/tiny_dataobj.p", 'rb')) #reload bc data object is mutable
    data_obj.workdir = absworkdir
    scraper2 = PhyscraperScrape(data_obj, ids)
    assert len(scraper2.data.aln) == 5

    scraper2.read_blast(blast_dir="tests/data/tiny_test_example/blast_files")
    scraper2.config.seq_len_perc = 0.998 #Change seq match percentage

    assert len(scraper2.new_seqs) == 40
    assert len(scraper2.new_seqs_otu_id) == 0

    scraper2.remove_identical_seqs()

    assert len(scraper2.new_seqs_otu_id) == 35
    sys.stdout.write("\n\nTest `remove_identical_seqs' passed\n\n")
except:
	sys.stdout.write("\n\nTest `remove_identical_seqs' FAILED\n\n")


#TODO - Check that out_dicts are correct,
#NO taxa labelled as not added, should be in the tree!!!
