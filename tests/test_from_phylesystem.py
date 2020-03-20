import physcraper
import sys
import pytest
from pytest import mark
from dendropy import DnaCharacterMatrix


def test_generate_ATT_from_phylesystem():
    study_id = "ot_350"
    tree_id = "Tr53296"
    seqaln = "tests/data/minitest.fas"
    mattype = "fasta"
    workdir = "tests/output/opentree"
    configfi = "tests/data/remotencbi.config"


    conf = physcraper.ConfigObj(configfi, interactive=False)
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

    data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                    workdir=workdir,
                                                    config_obj=conf,
                                                    study_id=study_id,
                                                    tree_id=tree_id)

    data_obj == True
    assert len(data_obj.tre.leaf_nodes())==16

test_generate_ATT_from_phylesystem()


@pytest.mark.xfail
def test_generate_ATT_from_phylesystem_fail():
    seqaln = "tests/data/input.fas"
    study_id = "pg_873"
    tree_id = "tree1679"
    seqaln = "tests/data/minitest.fas"
    mattype = "fasta"
    workdir = "tests/output/opentree"
    configfi = "tests/data/remotencbi.config"


    sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")

    conf = physcraper.ConfigObj(configfi, interactive=False)
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                    workdir=workdir,
                                                    config_obj=conf,
                                                    study_id=study_id,
                                                    tree_id=tree_id)
    