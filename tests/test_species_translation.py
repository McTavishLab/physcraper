import sys

from physcraper import get_mrca_ott, get_ott_taxon_info, OtuJsonDict, ConfigObj, IdDicts
from peyotl.sugar import taxomachine, tree_of_life

expected_json = {'otuSdoronicum': {'^ot:ottTaxonName': u'Senecio doronicum', '^physcraper:status': 'original', '^ot:ottId': 318436, '^user:TaxonName': 'Senecio_doronicum', '^ncbi:taxon': u'462523', '^ot:originalLabel': 'S_doronicum', '^physcraper:last_blasted': '1900/01/01'}, 'otuSlagascanus': {'^ot:ottTaxonName': u'Senecio lagascanus', '^physcraper:status': 'original', '^ot:ottId': 640718, '^user:TaxonName': 'Senecio_lagascanus', '^ncbi:taxon': u'1268580', '^ot:originalLabel': 'S_lagascanus', '^physcraper:last_blasted': '1900/01/01'}, 'otu2029doronicum': {'^ot:ottTaxonName': u'Senecio doronicum', '^physcraper:status': 'original', '^ot:ottId': 318436, '^user:TaxonName': 'Senecio_doronicum', '^ncbi:taxon': u'462523', '^ot:originalLabel': '2029_doronicum', '^physcraper:last_blasted': '1900/01/01'}, 'otuSlopezii': {'^ot:ottTaxonName': u'Senecio lopezii', '^physcraper:status': 'original', '^ot:ottId': 688688, '^user:TaxonName': 'Senecio_lopezii', '^ncbi:taxon': u'1268581', '^ot:originalLabel': 'S_lopezii', '^physcraper:last_blasted': '1900/01/01'}, 'otuSscopolii': {'^ot:ottTaxonName': u'Senecio scopolii', '^physcraper:status': 'original', '^ot:ottId': 688671, '^user:TaxonName': 'Senecio_scopolii', '^ncbi:taxon': u'1268589', '^ot:originalLabel': 'S_scopolii', '^physcraper:last_blasted': '1900/01/01'}}

sys.stdout.write("test_species_translation.py \n\n")

def test_species_translation():
	spn = "Mephitis mephitis"
	info = get_ott_taxon_info(spn)
	if info:
	    ottid, ottname, ncbi_id = info
	a = ottid == 231602

	tree_of_life.mrca(ott_ids=[ottid], wrap_response=False)


	ott_ids = [770315, 158484]
	ott_mrca = get_mrca_ott(ott_ids)
	b =  ott_mrca == 312031



	workdir="tests/output/tmp"
	configfi = "tests/data/test.config"
	id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
	otu_jsonfi = "{}/otu_dict.json".format(workdir)

	"""Tests if your own input files will generate a data object of class AlignTreeTax
	"""

	conf = ConfigObj(configfi, interactive=False)
	ids = IdDicts(conf, workdir=workdir)

	otu_json = OtuJsonDict(id_to_spn, ids)



	c = otu_json == expected_json
	assert a*b*c == 1