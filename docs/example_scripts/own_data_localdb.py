from physcraper import wrappers

seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
blacklist = None
workdir="tests/output/addLocal"
id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"

configfi = "tests/data/localblast.config"
threshold=10
selectby="blast" 
downto= None
ingroup_mrca = None
add_unpubl_seq = "tests/data/local_seqs"
id_to_spn_addseq = "tests/data/tipnTOspn_localAdd.csv"

wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         id_to_spn,
                         configfi,
                         selectby=selectby, 
                         downtorank=downto,
      			     ingroup_mrca=ingroup_mrca,
                         add_unpubl_seq=add_unpubl_seq,
                         id_to_spn_addseq_json=id_to_spn_addseq)
