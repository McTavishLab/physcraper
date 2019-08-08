import physcraper
from peyotl.sugar import tree_of_life, taxonomy, taxomachine

spp_name = "homo sapiens"
physcraper.opentree_helpers.get_ott_taxon_info("homo sapiens")

ottids = 515698,590452,643717
mrca = physcraper.get_mrca_ott(['515698','590452','643717'])

assert mrca == 1042120

resp = tree_of_life.mrca(ott_ids=['515698'], wrap_response=False)

expected_resp = {u'mrca': {u'taxon': {u'unique_name': u'Barnadesia', u'tax_sources': [u'ncbi:4237', u'gbif:8141286', u'irmng:1046392'], u'name': u'Barnadesia', u'rank': u'genus', u'ott_id': 515698}, u'supported_by': {u'pg_99@tree5885': u'_node1001466 ', u'ott3.0draft6': u'ott515698'}, u'was_constrained': True, u'node_id': u'ott515698', u'was_uncontested': True, u'resolves': {u'pg_2539@tree6294': u'_node1094881 ', u'ot_502@tree1': u'_node910 ', u'pg_1944@tree3959': u'_node740947 '}, u'num_tips': 14}, u'synth_id': u'opentree10.4', u'source_id_map': {u'pg_2539@tree6294': {u'git_sha': u'17ef8789bb8c7de13f285d5198c4e6c120695f80', u'study_id': u'pg_2539', u'tree_id': u'tree6294'}, u'pg_99@tree5885': {u'git_sha': u'17ef8789bb8c7de13f285d5198c4e6c120695f80', u'study_id': u'pg_99', u'tree_id': u'tree5885'}, u'ot_502@tree1': {u'git_sha': u'17ef8789bb8c7de13f285d5198c4e6c120695f80', u'study_id': u'ot_502', u'tree_id': u'tree1'}, u'ott3.0draft6': {u'taxonomy': u'ott3.0draft6'}, u'pg_1944@tree3959': {u'git_sha': u'17ef8789bb8c7de13f285d5198c4e6c120695f80', u'study_id': u'pg_1944', u'tree_id': u'tree3959'}}}

assert resp == expected_resp