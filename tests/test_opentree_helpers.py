from physcraper import treetaxon, opentree_helpers
import sys
import pytest
import xml

def test_root_tree_from_synth():
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='ott')
    opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='synth')

def test_get_dataset_from_treebase():
    # tb_id = 11681  # the dataset from this treebase id was giving trouble
    # treebase study S11681 corresponds to OpenTree study id pg_1976
    # url = "https://raw.githubusercontent.com/TreeBASE/supertreebase/master/data/treebase/S{}.xml".format(tb_id)
    study_id = "pg_1976"
    opentree_helpers.get_dataset_from_treebase(study_id)
