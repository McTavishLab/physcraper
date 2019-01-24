from physcraper import wrappers

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir="docs/example_scripts/output/OToL_filter"
configfi = "tests/data/localblast.config"


threshold = 2  # amount of sequences being kept by FilterBlast
selectby = "blast"  # how to select sequences in FilterBlast, either "length" or "blast"

ingroup_mrca = None  # must be OToL ID
shared_blast_folder = None # location to share blast runs across runs, see documentation

downtorank = None  # define filter rank, e.g. "species", "genus", if not defined, goes down to var/subsp
blacklist = None  # list with accession numbers, e.g. [XXX.1, YYY.1]
add_unpubl_seq = None
id_to_spn_addseq_json = None


## function to filter the blast results, if you want to keep all sequences found by blast, use standard_run()
wrappers.filter_OTOL(study_id,
                tree_id,
                seqaln,
		mattype,
                workdir,
                configfi,
                threshold,
                selectby=selectby,
                downtorank=downtorank,
                blacklist=blacklist,
                add_unpubl_seq=add_unpubl_seq,
                id_to_spn_addseq_json=id_to_spn_addseq_json,
                ingroup_mrca=ingroup_mrca,
                shared_blast_folder=shared_blast_folder)
