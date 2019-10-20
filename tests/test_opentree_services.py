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
    match = opentree_helpers.get_ottid_from_gbifid(gbif_id)
    assert match == 124230
    gbif_id2 = '66666'
    match2 = opentree_helpers.get_ottid_from_gbifid(gbif_id2)
    assert match2 == None

