
from physcraper import wrappers

seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "docs/example_scripts/output/own_standard_local"
configfi = "tests/data/localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"

ingroup_mrca = None
shared_blast_folder = None


# this function will keep all sequences found by blast which belong to the mrca,
# if you want to filter use filter_data_run()
wrappers.own_data_run(seqaln,
                  mattype,
                  trfn,
                  schema_trf,
                  workdir,
                  id_to_spn,
                  configfi,
                  ingroup_mrca=ingroup_mrca,
                  shared_blast_folder=shared_blast_folder)


