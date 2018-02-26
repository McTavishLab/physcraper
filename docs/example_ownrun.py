from physcraper import wrappers_nonstandard
#
seqaln= "docs/owndata/senecio_its.fasta"
mattype="fasta"
trfn= "docs/owndata/its_new.tre"
schema_trf = "newick"
workdir="docs/owndata/example_owndata_output"
otujson = "docs/owndata/ott_info_owntree.txt"
configfi = "example.config"
idtospname = "docs/owndata/uniquetip_to_name_its.csv"



wrappers_nonstandard.own_data_run(idtospname,
				 seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 otujson,
                 configfi)
