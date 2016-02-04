#!/usr/bin/env python
from physcraper import StudyInfo, PhyscraperSetup, PhyscraperScrape
import pickle


#LSU ASC tree example

study_id = "pg_238"
tree_id = "tree110"
seqaln = "rpb2.fas"
mattype="fasta"
runname="rpb2_test"


info = StudyInfo(study_id, tree_id, seqaln, mattype)

test = PhyscraperSetup(info, runname)
test.setup_physcraper()

pickle.dump(test, open('{}/{}_setup.p'.format(runname,runname), 'wb'))

test2 = PhyscraperScrape(test)
#test2.set_date("2016-01-13")
#test2 = pickle.load(open('{}/{}_scrape.p'.format(runname,runname),'rb'))
#test2.today = "2016-01-13"
test2.generate_streamed_alignment()
#test2.set_date("2016-01-15")
#test2.generate_streamed_alignment()
#test2.mrca_ncbi
#test2.run_blast()
#test2._blast_complete = 1
#test2.read_blast()
#test2.remove_identical_seqs()
#test2.write_query_seqs()
#test2.scrape()

#TODO need actual mini data set for test, need to fix dendropy pickle,
pickle.dump(test2, open('{}/{}_scrape.p'.format(runname,runname), 'wb'))
