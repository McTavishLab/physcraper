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
    configfi = "tests/data/test.config"


    data_obj = physcraper.generate_ATT_from_phylesystem(alnfile=seqaln,
                                                        aln_schema = mattype,
                                                        workdir=workdir,
                                                        configfile=configfi,
                                                        study_id=study_id,
                                                        tree_id=tree_id)

    data_obj == True
    assert(len(data_obj.tre.leaf_nodes())==16)
    assert(data_obj.otu_dict['Tl805422']['^ot:ottId'] ==  517518)
    assert(data_obj.otu_dict['Tl805431']['^ot:originalLabel']=='Dinemasporium_pseudostrigosum_CBS_717.85')
test_generate_ATT_from_phylesystem()


@pytest.mark.xfail
def test_generate_ATT_from_phylesystem_fail():
    seqaln = "tests/data/input.fas"
    study_id = "pg_873"
    tree_id = "tree1679"
    seqaln = "tests/data/minitest.fas"
    mattype = "fasta"
    workdir = "tests/output/opentree"
    configfi = "tests/data/test.config"


    sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")

    data_obj = physcraper.generate_ATT_from_phylesystem(alnfile=seqaln,
                                                        aln_schema = mattype,
                                                        workdir=workdir,
                                                        configfile=configfi,
                                                        study_id=study_id,
                                                        tree_id=tree_id)
    