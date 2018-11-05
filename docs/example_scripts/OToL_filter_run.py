from physcraper import wrappers

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
workdir="example_output"
configfi = "tests/data/test.config"

threshold = 2
selectby = "blast"
downtorank = "species"
ingroup_mrca = None

blacklist = None
add_unpubl_seq = None
id_to_spn_addseq_json = None
shared_blast_folder = None

wrappers.filter_OTOL(study_id,
                tree_id,
                seqaln,
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
