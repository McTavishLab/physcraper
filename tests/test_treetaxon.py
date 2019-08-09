import os
import filecmp
from physcraper import opentree_helpers
from physcraper.treetaxon import TreeTax


json_file = "tests/data/treetaxon/main.json"
exp_otu_dict = {u'Otuname11': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Phrynops', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 66456, '^ot:originalLabel': u'phrynops'}, u'Otuname10': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Emys orbicularis', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 733093, '^ot:originalLabel': u'emys_orbicularis'}, u'Otuname13': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Caretta', '^irmng:taxon': u'irmng:', '^worms:taxon': u'worms:', '^gbif:taxon': u'gbif:', '^ot:ottId': 66463, '^ot:originalLabel': u'caretta'}, u'Otuname12': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Caiman', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 335589, '^ot:originalLabel': u'caiman'}, u'Otuname15': {'^ncbi:taxon': u'ncbi:', '^ot:ottId': 284917, '^ot:originalLabel': u'chelonoidis_nigra', '^gbif:taxon': u'gbif:', '^ot:ottTaxonName': u'Chelonoidis nigra'}, u'Otuname14': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Python', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 675102, '^ot:originalLabel': u'python'}, u'Otuname16': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Podarcis', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 937560, '^ot:originalLabel': u'podarcis'}, u'Otuname9': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Alligator', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 335593, '^ot:originalLabel': u'alligator'}, u'Otuname8': {'^ncbi:taxon': u'ncbi:', '^irmng:taxon': u'irmng:', '^ot:ottId': 465090, '^ot:originalLabel': u'Xenopus', '^ot:ottTaxonName': u'Xenopus (genus in Deuterostomia)'}, u'Otuname1': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Protopterus', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 199350, '^ot:originalLabel': u'protopterus'}, u'Otuname3': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Gallus', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 153562, '^ot:originalLabel': u'Gallus'}, u'Otuname2': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Anolis', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 705358, '^ot:originalLabel': u'Anolis'}, u'Otuname5': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Monodelphis', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 122359, '^ot:originalLabel': u'Monodelphis'}, u'Otuname4': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Homo', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 770309, '^ot:originalLabel': u'Homo'}, u'Otuname7': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Taeniopygia', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 708325, '^ot:originalLabel': u'Taeniopygia'}, u'Otuname6': {'^ncbi:taxon': u'ncbi:', '^ot:ottTaxonName': u'Ornithorhynchus', '^irmng:taxon': u'irmng:', '^gbif:taxon': u'gbif:', '^ot:ottId': 962391, '^ot:originalLabel': u'Ornithorhynchus'}}
treefile = "tests/data/treetaxon/turtle.fa.1.treefile"
treeout = "tests/data/tmp/ottid.tre"
expected_tree = "tests/data/treetaxon/ottid.tre"


def test_load_bulk():
    otu_dict = opentree_helpers.bulk_tnrs_load(json_file)
    assert otu_dict == exp_otu_dict

def test_tree_taxon():
    tt = TreeTax(otu_json="tests/data/treetaxon/main.json", treefrom=treefile )
    tt.write_labelled(treeout, label = "^ot:ottId")
    assert os.path.isfile(treeout)
    assert filecmp.cmp(treeout, expected_tree)
    os.remove(treeout)

def test_get_synth_tree():
    ott_ids = set([exp_otu_dict[otu].get("^ot:ottId") for otu in exp_otu_dict])
    ott_ids = list(ott_ids)
    resp = opentree_helpers.get_tree_from_synth(ott_ids=ott_ids)
    print resp


def test_get_phyle_tree():
    tr = opentree_helpers.get_tree_from_study(study_id='pg_1144', tree_id='tree2324', label_format="name")


