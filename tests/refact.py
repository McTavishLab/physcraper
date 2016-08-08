from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts,  PhyscraperScrape
import dendropy
import pickle
import sys
import os


study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"



configfi = "tests/local.config"

conf = ConfigObj(configfi)

aln = dendropy.DnaCharacterMatrix.get(file=open(seqaln), schema=mattype)

data_obj = generate_ATT_from_phylesystem(aln,
                     "tmp",
                     study_id = study_id,
                     tree_id = tree_id,
                     phylesystem_loc = conf.phylesystem_loc)


data_obj.prune_short()
data_obj.write_files()
data_obj.write_labelled()


ids = IdDicts(conf, "tmp")

scraper =  PhyscraperScrape(data_obj, ids, conf)
scraper.run_blast()
scraper.read_blast()
scraper.remove_identical_seqs()
scraper.generate_streamed_alignment()

scraper.run_blast()
scraper.read_blast()
scraper.remove_identical_seqs()
scraper.generate_streamed_alignment()


#otu_json = "tests/minitest_otu.json"
#treefile = "tests/minitest.tre"
#info_obj2 = StudyInfo(seqaln,
#                     mattype,
                      #"tmp"
#                     otu_json = otu_json,
#                     treefile = treefile)

#data_obj2 = info_obj.generate_ATT()
'''
print("ROUND TWO")
data_obj2 = generate_ATT_from_files("tmp/physcraper.fas",
                                   mattype,
                                   "tmp",
                                   "tmp/physcraper.tre",
                                    "tmp/otus.json",
                                   ingroup_mrca=data_obj.ott_mrca)


scraper2 =  PhyscraperScrape(data_obj2, ids, conf)
scraper2.run_blast()
scraper2.read_blast()
'''