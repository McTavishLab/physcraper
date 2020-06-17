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
parser.add_argument("-d","--results_dir", help="results directory from run")
parser.add_argument("-t1","--original_tree", help="Original tree")
parser.add_argument("-t2","--updated_tree", help="Updated Tree")
parser.add_argument("-otu", "--otu_info", help="File with taxon information JSON")

schema = "newick"



args = parser.parse_args()

tns = dendropy.TaxonNamespace()


if args.results_dir:
    workdir = args.results_dir
    files = [f for f in os.listdir(workdir)]
    for file in files:
        if file.startswith('inputs_'):
            tag = file.split('.')[0].replace('inputs_', '')
    rundir = "{}/run_{}".format(workdir, tag)
    outputsdir = "{}/outputs_{}".format(workdir, tag)
    inputsdir = "{}/inputs_{}".format(workdir, tag)
    tree1_path = "{}/physcraper_{}.tre".format(inputsdir, tag)
    otu_json_path = "{}/otu_info_{}.json".format(rundir, tag)
    tree2_path = "{}/physcraper_{}.tre".format(outputsdir, tag)
else:
    tree1_path = args.original_tree
    tree2_path =- args.updated_tree
    otu_json_path = args.otu_info


tree1 = dendropy.Tree.get_from_path(tree1_path,
                                    schema,
                                    taxon_namespace=tns,
                                    preserve_underscores=True)
tree2 = dendropy.Tree.get_from_path(tree2_path,
                                    schema,
                                    taxon_namespace=tns,
                                    preserve_underscores=True)

otu_dict = json.load(open(otu_json_path, "r"))


leaves_t1 = set([leaf.taxon for leaf in tree1.leaf_nodes()])
leaves_t2 = set([leaf.taxon for leaf in tree2.leaf_nodes()])

new_tips = len(leaves_t1) - len(leaves_t2)
sys.stdout.write("{} new tips were added\n".format(new_tips))


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
sys.stdout.write("The RobinsonFoulds distance between these trees is {} and the weighted RF is {}".format(RF, weightedrf))


##Write out t2 for conflict with opentree
workdir = os.getcwd()


def conflict_tree(inputtree, otu_dict):
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
        return tmp_tree

tree_updated = conflict_tree(unpruned_tree2, otu_dict)
tree_orig = conflict_tree(tree1, otu_dict)


resp_updated = OT.conflict_str(tree_updated.as_string(schema="newick"), 'ott')
resp_orig = OT.conflict_str(tree_orig.as_string(schema="newick"), 'ott')
conflict = resp_updated.response_dict

for node in tree_updated:
    if node.taxon:
        node_id = node.taxon.label.split('_')[1]
        conf_node = conflict.get(node_id, False)
        if conf_node:
            new_label = "{} {}".format(conf_node['status'], conf_node['witness_name'])
            node.taxon.label = new_label
    else:
        node_id = node.label.split('_')[1]
        conf_node = conflict.get(node_id, False)
        if conf_node:
            new_label = "{} {}".format(conf_node['status'], conf_node['witness_name'])
            node.label = new_label

tree_updated.write(path = "conflict_label.tre", schema="newick")


'''
tree_updated_synth = conflict_tree(unpruned_tree2, otu_dict)
resp_updated_synth = OT.conflict_str(tree_updated_synth.as_string(schema="newick"), 'synth')
conflict_synth = resp_updated_synth.response_dict


for node in tree_updated_synth:
    if node.taxon:
        node_id = node.taxon.label.split('_')[1]
        conf_node = conflict_synth.get(node_id, False)
        if conf_node:
            new_label = "{} {}".format(conf_node['status'], conf_node['witness'])
            node.taxon.label = new_label
    else:
        node_id = node.label.split('_')[1]
        conf_node = conflict_synth.get(node_id, False)
        if conf_node:
            new_label = "{} {}".format(conf_node['status'], conf_node['witness'])
            node.label = new_label

'''