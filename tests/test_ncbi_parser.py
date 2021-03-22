from physcraper import ncbi_data_parser

from pytest import mark

slow = mark.slow


def test_ncbi_parser():
    ncbitax = ncbi_data_parser.Parser("./taxonomy/names.dmp","./taxonomy/nodes.dmp")

    taxid = ncbitax.get_id_from_name("Crassocephalum vitellinum")


    rankidgenus = ncbitax.get_downtorank_id(taxid, downtorank="genus")

    mrcaid = ncbitax.match_id_to_mrca(taxid, rankidgenus)

    mrca_id2 = ncbitax.match_id_to_mrca(17043521, taxid)
    
    name = ncbitax.get_name_from_id(1892268)

    rankid = ncbitax.get_downtorank_id(taxid)

    rank = ncbitax.get_rank(taxid)

    synonym = ncbitax.get_id_from_synonym("Elaps heterozonus")

    assert taxid == 1892268
    assert mrcaid == True
    assert mrca_id2 == False
    assert rankid == 1892268
    assert rankidgenus == 189210
    assert name == "Crassocephalum_vitellinum"
    assert rank == "species"
    assert synonym ==1117135



test_ncbi_parser()
	