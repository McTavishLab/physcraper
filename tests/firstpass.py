from physcraper import physcraper_setup, physcraper_scrape
import pickle


#LSU ASC tree example

study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/tree1679.fas"
mattype="fasta"
runname="asc_test_full"


test = physcraper_setup(study_id, tree_id, seqaln, mattype, runname)
test._read_config()
test._get_mrca()
assert(test.mrca_ott == 921280) #changes in the tree could actually change this...

test._reconcile_names()
test._prune()
test._write_files()
pickle.dump(test, open('{}/{}_setup.p'.format(runname,runname), 'wb'))

test2 = physcraper_scrape('{}/{}_setup.p'.format(runname,runname))
#test2.today = "2016-01-08"

#test2.mrca_ncbi
test2.run_blast()
#test2._blast_complete = 1
test2.read_blast()
test2.remove_identical_seqs()
test2.write_query_seqs()
#test2.scrape()
#TODO need actual mini data set for test, need to fix dendropy pickle,
pickle.dump(test2, open('{}/{}_scrape.p'.format(runname,runname), 'wb'))


#test3 = physcraper_add('{}_scrape.p'.format(test.runname))
