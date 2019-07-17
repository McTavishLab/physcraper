from physcraper import ncbi_data_parser

from pytest import mark

slow = mark.slow


def test_ncbi_parser():
    ncbitax = ncbi_data_parser.Parser("./tests/data/names.dmp","./tests/data/nodes.dmp")
    ncbitax.initialize()

    taxid = ncbitax.get_id_from_name("Crassocephalum vitellinum")


    rankidgenus = ncbitax.get_downtorank_id(taxid, downtorank="genus")

    mrcaid = ncbitax.match_id_to_mrca(taxid, rankidgenus)

    rankid = ncbitax.get_downtorank_id(taxid)

    rank = ncbitax.get_rank(taxid)

    synonym = ncbitax.get_id_from_synonym("Elaps heterozonus")

    assert taxid == 1892268
    assert mrcaid == True
    assert rankid == 1892268
    assert rankidgenus == 189210

    assert rank == "species"
    assert synonym ==1117135



test_ncbi_parser()
	