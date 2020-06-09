import sys
import os
import json
import argparse
import dendropy
import copy
import physcraper
from opentree import OT
from dendropy.calculate import treecompare
### Example
# python ../physcraper/bin/rf_distance.py -t1 pg_238/inputs_pg_238tree109_RPB2/physcraper_pg_238tree109_RPB2.tre -t2 pg_238/outputs_pg_238tree109_RPB2/physcraper_pg_238tree109_RPB2.tre -otu pg_238/run_pg_238tree109_RPB2/otu_info_pg_238tree109_RPB2.json 



parser = argparse.ArgumentParser()
parser.add_argument("-t1","--original_tree", help="Original tree")
parser.add_argument("-t2","--updated_tree", help="Updated Tree")
parser.add_argument("-otu", "--otu_info", help="File with taxon information JSON")

schema = "newick"



args = parser.parse_args()

tns = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get_from_path(args.original_tree,
                                    schema,
                                    taxon_namespace=tns,
                                    preserve_underscores=True)
tree2 = dendropy.Tree.get_from_path(args.updated_tree,
                                    schema,
                                    taxon_namespace=tns,
                                    preserve_underscores=True)

leaves_t1 = set([leaf.taxon for leaf in tree1.leaf_nodes()])
leaves_t2 = set([leaf.taxon for leaf in tree2.leaf_nodes()])

new_tips = len(leaves_t1) - len(leaves_t2)
sys.stdout.write("{} new tips were added\n".format(new_tips))

otu_dict = json.load(open(args.otu_info, "r"))

old_spp = set()
new_spp = set()

## Prune trees to same leaf set

for leaf in leaves_t1:
    species = otu_dict[leaf.label]['^ot:ottId']
    old_spp.add(species)

for leaf in leaves_t2:
    species = otu_dict[leaf.label]['^ot:ottId']
    new_spp.add(species)


sys.stdout.write("There were {} taxa in tree1, {} taxa in tree2\n".format(len(old_spp), len(new_spp)))


prune = leaves_t2.symmetric_difference(leaves_t1)

unpruned_tree2 = copy.deepcopy(tree2)
tree2.prune_taxa(prune)


RF = dendropy.calculate.treecompare.unweighted_robinson_foulds_distance(tree1, tree2)

#Weighted RF
weightedrf = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tree1, tree2)

tree2.write(path = "pruned_updated.tre", schema="newick")
print("The RobinsonFoulds distance between these trees is {} and the weighted RF is {}".format(RF, weightedrf))


##Write out t2 for conflict with opentree
workdir = os.getcwd()


def write_conflict_tree(inputtree, otu_dict):
        tmp_tree = copy.deepcopy(inputtree)
        new_names = set()
        i = 1
        for node in tmp_tree:
            i+=1
            if node.taxon:
                otu = otu_dict[node.taxon.label]
                ottid = otu['^ot:ottId']
                new_label = "_nd{}_ott{}".format(i, ottid)
                node.taxon.label = new_label
            else:
                node.label = "_nd{}_".format(i)
        return tmp_tree.as_string(schema="newick")

treestr_updated = write_conflict_tree(unpruned_tree2, otu_dict)
treestr_orig = write_conflict_tree(tree1, otu_dict)


resp_updated = OT.conflict_str(treestr_updated, 'ott')
resp_orig = OT.conflict_str(treestr_orig, 'ott')

print(resp_orig.response_dict)


