import sys

from physcraper import  OtuJsonDict, ConfigObj, IdDicts
from physcraper.opentree_helpers import get_mrca_ott, get_ott_taxon_info
from opentree import OT

def test_species_translation():
	spn = "Mephitis mephitis"
	info = get_ott_taxon_info(spn)
	if info:
	    ottid, ottname, ncbi_id = info
	assert ottid == 231602
	resp= get_mrca_ott(ott_ids=[ottid])


	ott_ids = [770315, 158484]
	ott_mrca = get_mrca_ott(ott_ids)
	assert  ott_mrca == 312031


def test_compare_json():
	expected_json = {'2029_doronicum': {'^ncbi:taxon': 462523, '^ot:ottTaxonName': 'Senecio doronicum', '^ot:ottId': 318436, '^ot:originalLabel': '2029_doronicum', '^user:TaxonName': 'Senecio_doronicum', '^physcraper:status': 'original', '^physcraper:last_blasted': None, '^physcraper:TaxonName': 'Senecio doronicum', '^ncbi:TaxonName': 'Senecio doronicum'}, 'S_doronicum': {'^ncbi:taxon': 462523, '^ot:ottTaxonName': 'Senecio doronicum', '^ot:ottId': 318436, '^ot:originalLabel': 'S_doronicum', '^user:TaxonName': 'Senecio_doronicum', '^physcraper:status': 'original', '^physcraper:last_blasted': None, '^physcraper:TaxonName': 'Senecio doronicum', '^ncbi:TaxonName': 'Senecio doronicum'}, 'S_lagascanus': {'^ncbi:taxon': 1268580, '^ot:ottTaxonName': 'Senecio lagascanus', '^ot:ottId': 640718, '^ot:originalLabel': 'S_lagascanus', '^user:TaxonName': 'Senecio_lagascanus', '^physcraper:status': 'original', '^physcraper:last_blasted': None, '^physcraper:TaxonName': 'Senecio lagascanus', '^ncbi:TaxonName': 'Senecio lagascanus'}, 'S_lopezii': {'^ncbi:taxon': 1268581, '^ot:ottTaxonName': 'Senecio lopezii', '^ot:ottId': 688688, '^ot:originalLabel': 'S_lopezii', '^user:TaxonName': 'Senecio_lopezii', '^physcraper:status': 'original', '^physcraper:last_blasted': None, '^physcraper:TaxonName': 'Senecio lopezii', '^ncbi:TaxonName': 'Senecio lopezii'}, 'S_scopolii': {'^ncbi:taxon': 1268589, '^ot:ottTaxonName': 'Senecio scopolii', '^ot:ottId': 688671, '^ot:originalLabel': 'S_scopolii', '^user:TaxonName': 'Senecio_scopolii', '^physcraper:status': 'original', '^physcraper:last_blasted': None, '^physcraper:TaxonName': 'Senecio scopolii', '^ncbi:TaxonName': 'Senecio scopolii'}}



	workdir="tests/output/tmp"
	configfi = "tests/data/test.config"
	id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
	otu_jsonfi = "{}/otu_dict.json".format(workdir)


	conf = ConfigObj(configfi)
	ids = IdDicts(configfi)

	otu_json = OtuJsonDict(id_to_spn, ids)

	print(otu_json)
	assert otu_json == expected_json

