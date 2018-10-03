from physcraper import wrappers

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"
workdir="example_output"
configfi = "tests/data/test.config"



wrappers.standard_run(study_id,
                      tree_id,
                      seqaln,
                      mattype,
                      workdir,
                      configfi)
