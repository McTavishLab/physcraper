import os
import filecmp
from physcraper import opentree_helpers
from physcraper import TreeTax


json_file = "tests/data/treetaxon/main.json"
exp_otu_dict = {u'Otuname11': {'^ncbi:taxon': u'8462', '^ot:ottTaxonName': u'Phrynops', '^irmng:taxon': u'1201383', '^gbif:taxon': u'2442114', '^physcraper:status': 'original', '^ot:ottId': 66456, '^ot:originalLabel': u'phrynops', '^physcraper:last_blasted': None}, u'Otuname10': {'^ncbi:taxon': u'82168', '^ot:ottTaxonName': u'Emys orbicularis', '^irmng:taxon': u'11010173', '^gbif:taxon': u'5220538', '^physcraper:status': 'original', '^ot:ottId': 733093, '^ot:originalLabel': u'emys_orbicularis', '^physcraper:last_blasted': None}, u'Otuname13': {'^ncbi:taxon': u'8466', '^ot:ottTaxonName': u'Caretta', '^irmng:taxon': u'1324374', '^worms:taxon': u'137066', '^gbif:taxon': u'2442177', '^physcraper:status': 'original', '^ot:ottId': 66463, '^ot:originalLabel': u'caretta', '^physcraper:last_blasted': None}, u'Otuname12': {'^ncbi:taxon': u'8497', '^ot:ottTaxonName': u'Caiman', '^irmng:taxon': u'1010136', '^gbif:taxon': u'5220195', '^physcraper:status': 'original', '^ot:ottId': 335589, '^ot:originalLabel': u'caiman', '^physcraper:last_blasted': None}, u'Otuname15': {'^ncbi:taxon': u'66189', '^ot:ottTaxonName': u'Chelonoidis nigra', '^gbif:taxon': u'5220266', '^physcraper:status': 'original', '^ot:ottId': 284917, '^ot:originalLabel': u'chelonoidis_nigra', '^physcraper:last_blasted': None}, u'Otuname14': {'^ncbi:taxon': u'37579', '^ot:ottTaxonName': u'Python', '^irmng:taxon': u'1031494', '^gbif:taxon': u'2454645', '^physcraper:status': 'original', '^ot:ottId': 675102, '^ot:originalLabel': u'python', '^physcraper:last_blasted': None}, u'Otuname16': {'^ncbi:taxon': u'42163', '^ot:ottTaxonName': u'Podarcis', '^irmng:taxon': u'1304163', '^gbif:taxon': u'2468993', '^physcraper:status': 'original', '^ot:ottId': 937560, '^ot:originalLabel': u'podarcis', '^physcraper:last_blasted': None}, u'Otuname9': {'^ncbi:taxon': u'8495', '^ot:ottTaxonName': u'Alligator', '^irmng:taxon': u'1039645', '^gbif:taxon': u'2441367', '^physcraper:status': 'original', '^ot:ottId': 335593, '^ot:originalLabel': u'alligator', '^physcraper:last_blasted': None}, u'Otuname8': {'^ncbi:taxon': u'8353', '^ot:ottTaxonName': u'Xenopus (genus in Deuterostomia)', '^irmng:taxon': u'1382944', '^physcraper:status': 'original', '^ot:ottId': 465090, '^ot:originalLabel': u'Xenopus', '^physcraper:last_blasted': None}, u'Otuname1': {'^ncbi:taxon': u'7885', '^ot:ottTaxonName': u'Protopterus', '^irmng:taxon': u'1295830', '^gbif:taxon': u'2441252', '^physcraper:status': 'original', '^ot:ottId': 199350, '^ot:originalLabel': u'protopterus', '^physcraper:last_blasted': None}, u'Otuname3': {'^ncbi:taxon': u'9030', '^ot:ottTaxonName': u'Gallus', '^irmng:taxon': u'1278118', '^gbif:taxon': u'2473720', '^physcraper:status': 'original', '^ot:ottId': 153562, '^ot:originalLabel': u'Gallus', '^physcraper:last_blasted': None}, u'Otuname2': {'^ncbi:taxon': u'28376', '^ot:ottTaxonName': u'Anolis', '^irmng:taxon': u'1301983', '^gbif:taxon': u'2468081', '^physcraper:status': 'original', '^ot:ottId': 705358, '^ot:originalLabel': u'Anolis', '^physcraper:last_blasted': None}, u'Otuname5': {'^ncbi:taxon': u'13615', '^ot:ottTaxonName': u'Monodelphis', '^irmng:taxon': u'1325350', '^gbif:taxon': u'7967492', '^physcraper:status': 'original', '^ot:ottId': 122359, '^ot:originalLabel': u'Monodelphis', '^physcraper:last_blasted': None}, u'Otuname4': {'^ncbi:taxon': u'9605', '^ot:ottTaxonName': u'Homo', '^irmng:taxon': u'1035772', '^gbif:taxon': u'2436435', '^physcraper:status': 'original', '^ot:ottId': 770309, '^ot:originalLabel': u'Homo', '^physcraper:last_blasted': None}, u'Otuname7': {'^ncbi:taxon': u'59728', '^ot:ottTaxonName': u'Taeniopygia', '^irmng:taxon': u'1265687', '^gbif:taxon': u'2493632', '^physcraper:status': 'original', '^ot:ottId': 708325, '^ot:originalLabel': u'Taeniopygia', '^physcraper:last_blasted': None}, u'Otuname6': {'^ncbi:taxon': u'9257', '^ot:ottTaxonName': u'Ornithorhynchus', '^irmng:taxon': u'1107086', '^gbif:taxon': u'2433375', '^physcraper:status': 'original', '^ot:ottId': 962391, '^ot:originalLabel': u'Ornithorhynchus', '^physcraper:last_blasted': None}}
treefile = "tests/data/treetaxon/turtle.fa.1.treefile"
treeout = "tests/data/tmp/ottid.tre"
expected_tree = "tests/data/treetaxon/ottid.tre"


def test_load_bulk():
    otu_dict = opentree_helpers.bulk_tnrs_load(json_file)
    assert otu_dict == exp_otu_dict

def test_tree_taxon():
    tt = TreeTax(otu_json="tests/data/treetaxon/main.json", treefrom=treefile )
    tt.write_labelled(path=treeout, label = "^ot:ottId", norepeats=False)
    assert os.path.isfile(treeout)
    assert filecmp.cmp(treeout, expected_tree)
    os.remove(treeout)

def test_get_synth_tree():
    ott_ids = set([exp_otu_dict[otu].get("^ot:ottId") for otu in exp_otu_dict])
    ott_ids = list(ott_ids)
    resp = opentree_helpers.get_tree_from_synth(ott_ids=ott_ids)
    print(resp)

def test_get_phyle_tree():
    tr = opentree_helpers.get_tree_from_study(study_id='pg_1144', tree_id='tree2324', label_format="name")

