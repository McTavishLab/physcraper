from physcraper import physcraper_setup, physcraper_scrape
import pickle


#LSU ASC tree example
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/tree1679.fas"
mattype="fasta"
runname="asc_test"


test = physcraper_setup(study_id, tree_id, seqaln, mattype, runname)
#test._read_config()
test._get_mrca()
assert(test.mrca_ott == 921280) #changes in the tree could actually change this...

test._reconcile_names()
test._prune()
test._write_files()
#pickle.dump(test, open('{}.p'.format(test.runname), 'wb'))

test2 = physcraper_scrape('{}_setup.p'.format(test.runname))
test2.mrca_ncbi
test2.run_blast()

#TODO need actual mini data set for test, need to fix dendropy pickle,