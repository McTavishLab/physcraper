from physcraper import wrappers

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "docs/example_scripts/output/OToL_standard"
configfi = "tests/data/localblast.config"

wrappers.standard_run(study_id,
                      tree_id,
                      seqaln,
                      mattype,
                      workdir,
                      configfi)
