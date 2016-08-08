from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts,  PhyscraperScrape
from dendropy import DnaCharacterMatrix
import pickle
import sys
import os

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"

#
seqaln = "tests/data/minitest.fas"
mattype="fasta"

workdir="example"

#A configuration file stores the key config parameters, such as teh location of phylesystem, and BLAST
configfi = "example.config"
#read the config file into a configuration object
conf = ConfigObj(configfi)

aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

#Generate an linked Alignment-Tree-Taxa object
data_obj = generate_ATT_from_phylesystem(aln=aln,
                     workdir=workdir,
                     study_id = study_id,
                     tree_id = tree_id,
                     phylesystem_loc = conf.phylesystem_loc)



#Alternately, the alignemnt tacon tree object can be generated from a local tree using "generate_ATT_from_files"
#This requires taxon label amppings for teh tip labels on the tree.

#Prune sequnces below a certain length threshold
#This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
data_obj.prune_short()

data_obj.write_files()
data_obj.write_labelled()


#Mapping identifiers between OpenTree and NCBI requires and identifier dict object
ids = IdDicts(conf, workdir=workdir)


#Now combine the data, the ids, and the configuration into a single physcraper scrape object
scraper =  PhyscraperScrape(data_obj, ids, conf)

#run the ananlyses
scraper.run_blast()
scraper.read_blast()
scraper.remove_identical_seqs()
scraper.generate_streamed_alignment()


#Keep running it as long as you are pulling down new sequences
while len(scraper.new_seqs) > 0: 
  scraper.run_blast()
  scraper.read_blast()
  scraper.remove_identical_seqs()
  scraper.generate_streamed_alignment()

