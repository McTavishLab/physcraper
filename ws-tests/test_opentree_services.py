import physcraper
from peyotl.sugar import tree_of_life, taxonomy, taxomachine

spp_name = "homo sapiens"
physcraper.get_ott_taxon_info("homo sapiens")

ottids = 515698,590452,643717
mrca = physcraper.get_mrca_ott(['515698','590452','643717'])

assert mrca == 46248
tree_of_life.mrca(ott_ids=['515698'], wrap_response=False)