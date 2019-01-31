import sys

from physcraper import get_mrca_ott, get_ott_taxon_info, OtuJsonDict, ConfigObj, IdDicts
from peyotl.sugar import taxomachine, tree_of_life


def test_species_translation():
	spn = "Mephitis mephitis"
	info = get_ott_taxon_info(spn)
	if info:
	    ottid, ottname, ncbi_id = info
	assert ottid == 231602

	tree_of_life.mrca(ott_ids=[ottid], wrap_response=False)


	ott_ids = [770315, 158484]
	ott_mrca = get_mrca_ott(ott_ids)
	assert  ott_mrca == 312031


def test_compare_json():
	expected_json = {'otuSdoronicum':
						 {'^ncbi:taxon': u'462523', '^ot:ottTaxonName': u'Senecio doronicum', '^ncbi:TaxonName': 'Senecio doronicum', '^physcraper:TaxonName': 'Senecio doronicum', '^physcraper:status': 'original', '^ot:ottId': 318436, '^user:TaxonName': 'Senecio_doronicum', '^ot:originalLabel': 'S_doronicum', '^physcraper:last_blasted': '1900/01/01'}, 
					'otuSlagascanus': 
						{'^ncbi:taxon': u'1268580', '^ot:ottTaxonName': u'Senecio lagascanus', '^ncbi:TaxonName': 'Senecio lagascanus', '^physcraper:TaxonName': 'Senecio lagascanus', '^physcraper:status': 'original', '^ot:ottId': 640718, '^user:TaxonName': 'Senecio_lagascanus', '^ot:originalLabel': 'S_lagascanus', '^physcraper:last_blasted': '1900/01/01'}, 
					'otu2029doronicum': 
						{'^ncbi:taxon': u'462523', '^ot:ottTaxonName': u'Senecio doronicum', '^ncbi:TaxonName': 'Senecio doronicum', '^physcraper:TaxonName': 'Senecio doronicum', '^physcraper:status': 'original', '^ot:ottId': 318436, '^user:TaxonName': 'Senecio_doronicum', '^ot:originalLabel': '2029_doronicum', '^physcraper:last_blasted': '1900/01/01'}, 
					'otuSlopezii': 
						{'^ncbi:taxon': u'1268581', '^ot:ottTaxonName': u'Senecio lopezii', '^ncbi:TaxonName': 'Senecio lopezii', '^physcraper:TaxonName': 'Senecio lopezii', '^physcraper:status': 'original', '^ot:ottId': 688688, '^user:TaxonName': 'Senecio_lopezii', '^ot:originalLabel': 'S_lopezii', '^physcraper:last_blasted': '1900/01/01'}, 
					'otuSscopolii': 
						{'^ncbi:taxon': u'1268589', '^ot:ottTaxonName': u'Senecio scopolii', '^ncbi:TaxonName': 'Senecio scopolii', '^physcraper:TaxonName': 'Senecio scopolii', '^physcraper:status': 'original', '^ot:ottId': 688671, '^user:TaxonName': 'Senecio_scopolii', '^ot:originalLabel': 'S_scopolii', '^physcraper:last_blasted': '1900/01/01'}
					}


	workdir="tests/output/tmp"
	configfi = "tests/data/test.config"
	id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
	otu_jsonfi = "{}/otu_dict.json".format(workdir)



	conf = ConfigObj(configfi, interactive=False)
	ids = IdDicts(conf, workdir=workdir)

	otu_json = OtuJsonDict(id_to_spn, ids)

	print(otu_json)
	assert otu_json == expected_json
	
