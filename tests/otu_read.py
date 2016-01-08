from physcraper import physcraper_setup, physcraper_scrape
import pickle


#LSU ASC tree example

study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/tree1679.fas"
mattype="fasta"
runname="asc_test"


test = physcraper_setup(study_id, tree_id, seqaln, mattype, runname)
