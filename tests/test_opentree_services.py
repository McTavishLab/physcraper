import physcraper
from peyotl.sugar import tree_of_life, taxonomy, taxomachine
from pytest import mark
import pytest
from physcraper import opentree_helpers

def test_opentree_service():
    spp_name = "homo sapiens"
    opentree_helpers.get_ott_taxon_info("homo sapiens")
    ottids = 515698,590452,643717
    mrca = physcraper.get_mrca_ott(['515698','590452','643717'])
    assert mrca == 1042120
    resp = tree_of_life.mrca(ott_ids=['515698'], wrap_response=False)
    assert resp['mrca']['taxon']['unique_name'] == 'Barnadesia'
    assert opentree_helpers.check_if_ottid_in_synth(1) == 0
    assert opentree_helpers.check_if_ottid_in_synth(878252) == 0
    assert opentree_helpers.check_if_ottid_in_synth(983181) == 1
    opentree_helpers.check_if_ottid_in_synth("10000000000")


def test_ottids_from_gbifids():
    gbif_id = 2440447
    match = opentree_helpers.get_ottids_from_gbifids(gbif_id)
    assert len(match) == 1
    assert match[2440447] == 124230
    gbif_ids = ['1428475', '2440447', '2235372', '66666']
    matches = opentree_helpers.get_ottids_from_gbifids(gbif_ids)
    assert len(matches) == 4
    assert matches[2440447] == 124230
    assert matches[66666] == None

