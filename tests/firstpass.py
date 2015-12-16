from physcraper import physcraper_setup



#LSU ASC tree example
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tree1679.fas"
mattype="fasta"
runname="asc_test"


test = physcraper_setup(study_id, tree_id, seqaln, mattype, runname)
#test._read_config()
test._get_mrca()
assert(test.mrca == 921280) #changes in the tree could actually change this...


