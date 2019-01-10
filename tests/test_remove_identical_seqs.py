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

    print("start")
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.ids.otu_rank = {}
    scraper.config.gifilename = False

    scraper._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    #scraper.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    scraper.read_blast_wrapper(blast_dir=blast_dir)

    a = len(scraper.new_seqs) == 40
    b = len(scraper.data.aln) == 5
    c = len(scraper.new_seqs_otu_id) == 0

    scraper.remove_identical_seqs()

    d = len(scraper.new_seqs) == 40
    e = len(scraper.data.aln) == 5
    f = len(scraper.new_seqs_otu_id) == 38
    g = 1
    for taxon in scraper.data.tre.taxon_namespace:
        h = taxon.label in scraper.data.otu_dict
        g = g*h
        status = scraper.data.otu_dict[taxon.label].get(u'^physcraper:status')
        i = status in ('original', 'query')
        g = g*i

    # Second test checks that seq len prec is affecting results
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb')) #reload bc data object is mutable
    data_obj.workdir = absworkdir
    scraper2 = PhyscraperScrape(data_obj, ids)
    scraper2.ids.otu_rank = {}

    scraper2.config.gifilename = False
    j = len(scraper2.data.aln) == 5
    # scraper2.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    scraper2.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
    scraper2.config.seq_len_perc = 0.998  # Change seq len percentage from default of 75%

    k = len(scraper2.new_seqs) == 40
    l = len(scraper2.new_seqs_otu_id) == 0

    scraper2.remove_identical_seqs()
    # print(scraper2.data.otu_dict)
    # print(len(scraper.new_seqs_otu_id), 38)
    # print(len(scraper2.new_seqs_otu_id), 36)
    m = len(scraper2.new_seqs_otu_id) == 36
    count = 0
    assert a*b*c*d*e*f*g*h*i*j*k*l*m ==True
  