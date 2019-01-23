
from physcraper import wrappers

seqaln = "tests/data/tiny_test_example/test.fas"  # alignment file
mattype = "fasta" # format of alignment
trfn = "tests/data/tiny_test_example/test.tre"  # tree file
schema_trf = "newick"  # format of tree file
workdir = "docs/example_scripts/output/own_data_filter"  # working directory
configfi = "tests/data/localblast.config"  # path to config file
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"  # path to file where phylogeny tip names corresponding to taxon names

threshold = 2  # amount of sequences being kept by FilterBlast
selectby = "blast"  # how to select sequences in FilterBlast, either "length" or "blast"

ingroup_mrca = None  # must be OToL ID
shared_blast_folder = None # location to share blast runs across runs, see documentation
downtorank = None  # define filter rank, e.g. "species", "genus"
blacklist = None  # list with accession numbers, e.g. [XXX.1, YYY.1]
add_unpubl_seq = None  # path to folder with unpublished sequences in fasta format
id_to_spn_addseq_json = None # path to file where sequence names correspond to taxon names

# function to filter the blast results, 
# if you want to keep all sequences found by blast, use own_data_run()
wrappers.filter_data_run(seqaln,
                     mattype,
                     trfn,
                     schema_trf,
                     workdir,
                     threshold,
                     id_to_spn,
                     configfi,
                     downtorank=downtorank,
                     selectby=selectby,
                     blacklist=blacklist,
                     add_unpubl_seq=add_unpubl_seq,
                     id_to_spn_addseq_json=id_to_spn_addseq_json,
                     ingroup_mrca=ingroup_mrca,
                     shared_blast_folder=shared_blast_folder
                     )

