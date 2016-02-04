#!/usr/bin/env python
from physcraper import StudyInfo, PhyscraperSetup, PhyscraperScrape
import pickle
import sys

study_id = "pg_238"
tree_id = "tree110"
seqaln = "shard/ascomycota_rpb1/rpb1.fas"
mattype="fasta"
runname="rpb1"
 
sys.stdout.write("setting up StudyINfo\n")
sys.stdout.flush()

info = StudyInfo(study_id, tree_id, seqaln, mattype)
sys.stdout.write("setting up scrape instance\n")
sys.stdout.flush()
 
test = PhyscraperSetup(info, runname, configfi="config")
test.setup_physcraper()
sys.stdout.write("Instance set up\n")
sys.stdout.flush()

pickle.dump(test, open('{}/{}_setup.p'.format(runname,runname), 'wb'))
 
test2 = PhyscraperScrape(test)
test2.generate_streamed_alignment()

#TODO need actual mini data set for test, need to fix dendropy pickle,
pickle.dump(test2, open('{}/{}_scrape.p'.format(runname,runname), 'wb'))