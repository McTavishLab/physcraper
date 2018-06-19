import pickle
import sys
import os
from physcraper import ConfigObj, PhyscraperScrape, IdDicts

#Function we want to test is scrape.remove_identical_seqs()
#What are the inputs?  
#physracper.scrape object, with new sequences read in.
#to make that we need: input data, idObject, and a configuration object.

##todo Make Sure 
sys.stdout.write("Running test remove_identical_seqs\n\n")
workdir="tests/data/tmp/owndata"
absworkdir = os.path.abspath(workdir)
conf = ConfigObj("tests/data/test.config")



try:
   data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
   data_obj.workdir = absworkdir
   ids = IdDicts(conf, workdir=data_obj.workdir)
   ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
   sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
   sys.stdout.write("\n\nTest `remove_identical_seqs' FAILED\n\n")
   sys.exit()


scraper =  PhyscraperScrape(data_obj, ids)
scraper._blasted = 1
blast_dir = "tests/data/tiny_test_example/blast_files"
scraper.read_blast(blast_dir=blast_dir)

a = len(scraper.new_seqs) == 40
b = len(scraper.data.aln) == 5
c =  len(scraper.new_seqs_otu_id) == 0

scraper.remove_identical_seqs()

d = len(scraper.new_seqs) == 40
e = len(scraper.data.aln) == 5
f = len(scraper.new_seqs_otu_id) == 38


g = 1
for taxon in scraper.data.tre.taxon_namespace:
    h = taxon.label in scraper.data.otu_dict
    g = g*h
    status =  scraper.data.otu_dict[taxon.label].get(u'^physcraper:status')
    i = status in ('original', 'query')
    g = g*i

#Second test checks that seq len prec is affecting results
data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb')) #reload bc data object is mutable
data_obj.workdir = absworkdir
scraper2 = PhyscraperScrape(data_obj, ids)
j = len(scraper2.data.aln) == 5

scraper2.read_blast(blast_dir="tests/data/precooked/tiny_test_example/blast_files")
scraper2.config.seq_len_perc = 0.998 #Change seq match percentage

k = len(scraper2.new_seqs) == 40
l = len(scraper2.new_seqs_otu_id) == 0

scraper2.remove_identical_seqs()

m = len(scraper2.new_seqs_otu_id) == 37


if a*b*c*d*e*f*g*h*i*j*k*l*m:
    sys.stdout.write("\n\nTest `remove_identical_seqs' passed\n\n")
else:
    sys.stderr.write("\n\nTest `remove_identical_seqs' FAILED\n\n")