import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, debug
from dendropy import DnaCharacterMatrix
from copy import deepcopy#


#
workdir_its = "./tiny_comb_its"
workdir_ets = "./tiny_comb_ets"


pickle_fn = "scrape_checkpoint.p"


workdir_comb = "./tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, "ets": {"workdir": workdir_ets, "pickle": pickle_fn}}

##############
# print("{}/{}".format(workdir, pickle_fn))
# scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn),'rb'))
# print(scrape)


#############

# print(workdir_comb)

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
						user_concat=None)

print(type(conc))














#
#
# # print(some)
#
# if conc.short_concat_seq:
#     print(conc.short_concat_seq)
#     print(conc.concatenated_aln.taxon_namespace)
#     conc.concatenated_aln.remove_sequences(conc.short_concat_seq)
#     conc.tre_as_start.prune_taxa(conc.short_concat_seq)
#     conc.tre_as_start.prune_taxa_with_labels(conc.short_concat_seq)#sometimes it does not delete it with the statement before. Tried to figure out why, have no clue yet.
#     #conc.aln.taxon_namespace.remove_taxon_label(tax)
#     aln_ids = set()
#     for tax in conc.concatenated_aln:
#         aln_ids.add(tax.label)
#     for leaf in conc.tre_as_start.leaf_nodes():
#         if leaf.taxon.label not in aln_ids:
#             conc.tre_as_start.prune_taxa([leaf])
#             conc.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
#             conc.tre_as_start.prune_taxa_with_labels([leaf])
# # print(tre_as_start.as_string(schema="newick"))
# # print(conc.concatenated_aln.as_string(schema="fasta"))

### maybe len of taxonnamespace is not reduced because i did not delete the nodes?
print(len(conc.tre_as_start.taxon_namespace))
print(len(conc.concatenated_aln.taxon_namespace))
print(conc.tre_as_start.taxon_namespace)
print(conc.concatenated_aln.taxon_namespace)

tre_ids = set()
for tax in conc.tre_as_start.taxon_namespace:
	tre_ids.add(tax.label)


aln_ids = set()
for tax in conc.concatenated_aln.taxon_namespace:
	aln_ids.add(tax.label)

# debug(len(conc.otu_dict.keys()))
# debug(len(aln_ids))
debug([item for item in tre_ids if item not in aln_ids])
debug([item for item in aln_ids if item not in tre_ids])

# # print(some)
# print("writing files")
# tre_as_start_str = conc.tre_as_start.as_string(schema="newick",
#                             unquoted_underscores=True,
#                             suppress_rooting=True)
#
# fi = open("{}/{}".format(conc.workdir, "starting.tre"), "w")
# fi.write(tre_as_start_str)
# fi.close()
# conc.concatenated_aln.write(path="{}/{}".format(conc.workdir, "concat_red.fasta"), schema="fasta")
conc.place_new_seqs()
conc.est_full_tree()