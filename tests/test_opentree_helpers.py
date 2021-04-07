from physcraper import treetaxon, opentree_helpers
import sys
import pytest

def test_root_tree_from_synth():
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='ott')
    opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='synth')
