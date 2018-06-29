import sys
import os
import json
import pickle
import random
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, Concat
from dendropy import DnaCharacterMatrix
#


#
workdir_its = "tiny_comb_its"
workdir_ets = "tiny_comb_ets"


pickle_fn = "scrape_checkpoint.p"


workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, "ets": {"workdir": workdir_ets, "pickle": pickle_fn}}

##############
# print("{}/{}".format(workdir, pickle_fn))
# scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn),'rb'))
# print(scrape)


#############

# print(workdir_comb)

concat = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb, 
						user_concat=None)



# print(prune_shortest)


####i"m here
#shortest seq need to be pruned from aln and the starting tree for reconstruction
# first = True
# len_all_taxa = {}
# for gene in concat.single_runs:
# 	# if first:

# 	# 	len_aln = len(concat.single_runs[gene].data.aln.taxon_namespace)
# 	# 	len2 =0
# 	# 	first = False
# 	# else:
# 	# 	len2 = len(concat.single_runs[gene].data.aln.taxon_namespace)
# 	len_aln_taxa = len(concat.single_runs[gene].data.aln.taxon_namespace)
# 	len_all_taxa[gene] = len_aln_taxa

# for gene, len_item in len_all_taxa.items():
# 	if first:
# 		len_max = len_item
# 		gene_max = gene
# 		first = False
# 	if len_item > len_max:
# 		len_max = len_item
# 		gene_max = gene

# # print(gene_max, len_max)



# # rewrite otus of tre as start into spn:
# # for otu in concat.single_runs[gene_max].data.tre.taxon_namespace:
# #     data = concat.single_runs[gene_max].data.otu_dict[otu.label]
# #     spn = data['^ot:ottTaxonName']
# #     otu.label = spn






# tre_as_start = concat.single_runs[gene_max].data.tre
tre_len_before = len(concat.tre_as_start.taxon_namespace)
print(len(concat.tre_as_start.taxon_namespace))
print(len(concat.concatenated_aln.taxon_namespace))
#  from prune short:
if concat.short_concat_seq:
    debug(concat.short_concat_seq)
    print(concat.concatenated_aln.taxon_namespace)
    concat.concatenated_aln.remove_sequences(concat.short_concat_seq)
    tre_as_start.prune_taxa(concat.short_concat_seq)
    tre_as_start.prune_taxa_with_labels(concat.short_concat_seq)#sometimes it does not delete it with the statement before. Tried to figure out why, have no clue yet.
    #concat.aln.taxon_namespace.remove_taxon_label(tax)
    aln_ids = set()
    for tax in concat.concatenated_aln:
        aln_ids.add(tax.label)
    for leaf in tre_as_start.leaf_nodes():
            if leaf.taxon.label not in aln_ids:
                concat.tre_as_start.prune_taxa([leaf])
                concat.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                concat.tre_as_start.prune_taxa_with_labels([leaf])
# print(tre_as_start.as_string(schema="newick"))
# print(concat.concatenated_aln.as_string(schema="fasta"))

### maybe len of taxonnamespace is not reduced because i did not delete the nodes?
print(len(tre_as_start.taxon_namespace))
print(len(concat.concatenated_aln.taxon_namespace))
print(tre_as_start.taxon_namespace)
print(concat.concatenated_aln.taxon_namespace)
print(some)
tre_as_start_str = tre_as_start.as_string(schema="newick",
                            unquoted_underscores=True,
                            suppress_rooting=True)

fi = open("{}/{}".format(concat.workdir, "starting.tre"), "w")
fi.write(tre_as_start_str)
fi.close()
concat.concatenated_aln.write(path="{}/{}".format(concat.workdir, "concat_red.fasta"), schema="fasta")
concat.place_new_seqs()
concat.est_full_tree()