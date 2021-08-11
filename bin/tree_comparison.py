#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  Physcraper Library.
##
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Luna L. Sanchez Reyes, Martha Kandziora, Emily Jane McTavish. (2020).
##     Physcraper: A Python package for continual update of evolutionary
##     estimates using the Open Tree of Life. BioRxiv 2020.09.15.299156.
##     doi: doi.org/10.1101/2020.09.15.299156
##
##############################################################################

"""
This module handles the code to generate the Physcraper command line interface
to compare two trees using (1) the Open Tree of Life conflict functions, and
(2) using the Robinson-Foulds statistic calculated using Dendropy functions.
"""

import sys
import os
import json
# import subprocess
import argparse
import copy
import dendropy
from dendropy.calculate import treecompare
import physcraper
from physcraper.opentree_helpers import root_tree_from_synth, conflict_tree, ottids_in_synth
from opentree import OT
### Example
# python ../physcraper/bin/rf_distance.py
# -t1 pg_238/inputs_pg_238tree109_RPB2/physcraper_pg_238tree109_RPB2.tre
# -t2 pg_238/outputs_pg_238tree109_RPB2/physcraper_pg_238tree109_RPB2.tre
# -otu pg_238/run_pg_238tree109_RPB2/otu_info_pg_238tree109_RPB2.json



parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="Physcraper output directory from a single run")
parser.add_argument("-t1", "--original_tree", help="Original tree")
parser.add_argument("-t2", "--updated_tree", help="Updated Tree")
parser.add_argument("-otu", "--otu_info", help="File with taxon information in JSON format")
parser.add_argument("-og", "--outgroup", nargs='+', help="otu ids of outgroup taxa for rooting")
parser.add_argument("-o", "--outputdir", help="Name for output directory")


schema = "newick"

args = parser.parse_args()


try:
    assert args.outputdir
except AssertionError:
    sys.stderr.write("ERROR: Output directory (-o) is required.\n")
    sys.exit(-1)


comparisondir = args.outputdir
if not os.path.exists(comparisondir):
    os.mkdir(comparisondir)


if args.results_dir:
    assert(os.path.exists(args.results_dir)), "Results directory {} not found\n".format(args.results_dir)
else:
    assert(os.path.exists(args.original_tree)), "Original tree {} not found\n".format(args.original_tree)
    assert(os.path.exists(args.updated_tree)), "Updated tree {} not found\n".format(args.updated_tree)
    assert(os.path.exists(args.otu_info)), "Otu_info {} not found\n".format(args.otu_info)


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
    tree2_path = args.updated_tree
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


leaves_t1 = {leaf.taxon.label for leaf in tree1.leaf_nodes()}
leaves_t2 = {leaf.taxon.label for leaf in tree2.leaf_nodes()}

old_spp = set()
new_spp = set()


for leaf in leaves_t1:
    species = otu_dict[leaf]['^ot:ottId']
    old_spp.add(species)

for leaf in leaves_t2:
    species = otu_dict[leaf]['^ot:ottId']
    new_spp.add(species)

if None in old_spp:
    old_spp.remove(None)


if None in new_spp:
    new_spp.remove(None)

new_tips = len(leaves_t2) - len(leaves_t1)
sys.stdout.write("{} new tips were added\n".format(new_tips))
tree2.write(path = "{}/before_rooting.tre".format(comparisondir), schema="newick")
#['otu376420','otu376439','otu376452']
if args.outgroup:
    outgroup = args.outgroup
    sys.stdout.write("Rooting tree using {} as outgroup\n".format(", ".join(outgroup)))
    for tip in outgroup:
        assert(tip in leaves_t2), "label {} not found in updated tree {}\n".format(tip, tree2_path)
    mrca = tree2.mrca(taxon_labels=outgroup)
    tree2.reroot_at_node(mrca, update_bipartitions=True)
else:
    # try:
    rooted = root_tree_from_synth(tree2, otu_dict, base='ott')
    ## In which cases will auto-root fail????
    # except:
        # sys.stdout.write("Auto-rooting failed, conflict results may be spurious.\n")

##Write out tree2 for conflict with opentree
tree2.write(path = "{}/after_rooting.tre".format(comparisondir), schema="newick")



## THIS SECTION LOOKS AT TAXA

physcraper_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
synthfile = open("{}/taxonomy/ottids_in_synth.txt".format(physcraper_dir))
ottids_in_synth = ottids_in_synth()


sys.stdout.write("\nThere were {} new taxa in the updated tree\n".format(len(new_spp) - len(old_spp)))
msg1 = "Of the {} taxa in original tree {} are not included in synthesis phylogenies.\n"
sys.stdout.write(msg1.format(len(old_spp), len(old_spp.difference(ottids_in_synth))))
msg2 = "Of the {} taxa in updated tree {} are not included in synthesis phylogenies.\n\n"
sys.stdout.write(msg2.format(len(new_spp), len(new_spp.difference(ottids_in_synth))))


ids = physcraper.IdDicts()
sys.stdout.write("Taxa with only taxonomic information in the OpenTree synthetic tree (so far!) are:\n")
for tax in new_spp.difference(ottids_in_synth):
    taxname = ids.ott_to_name.get(tax, '-')
    sys.stdout.write("ott{}: {}\n".format(tax, taxname))

## This section does tree comparison
## Prune trees to same leaf set

prune = leaves_t2.symmetric_difference(leaves_t1)

unpruned_tree1 = copy.deepcopy(tree1)
unpruned_tree2 = copy.deepcopy(tree2)
tree2.prune_taxa(prune)


RF = treecompare.unweighted_robinson_foulds_distance(tree1, tree2)

#Weighted RF
weightedrf = treecompare.weighted_robinson_foulds_distance(tree1, tree2)


for tax in tns:
    if tax.label in otu_dict:
        tax.label = tax.label + "_" + str(otu_dict[tax.label].get('^ot:ottTaxonName'))

## write put with tip labels that have taxon names
tree1.write(path = "{}/original.tre".format(comparisondir), schema="newick")
tree2.write(path = "{}/pruned_updated.tre".format(comparisondir), schema="newick")
msg3 = "\n\nThe RF distance between the matched tips in the trees is {} and the weighted RF is {}\n"
sys.stdout.write(msg3.format(RF, weightedrf))

workdir = comparisondir
tree_updated = conflict_tree(unpruned_tree2, otu_dict)
tree_orig = conflict_tree(unpruned_tree1, otu_dict)


resp_orig = OT.conflict_str(tree_orig.as_string(schema="newick"), 'ott')
resp_updated = OT.conflict_str(tree_updated.as_string(schema="newick"), 'ott')


conflict_orig = resp_orig.response_dict


orig_conf_taxa = set()
for node in conflict_orig:
    if conflict_orig[node]['status'] == 'conflicts_with':
        witness = conflict_orig[node]['witness_name']
        for tax in witness:
            orig_conf_taxa.add(tax)

sys.stdout.write("\nOriginal tree conflicts with {} taxa in the OpenTree taxonomy:\n".format(len(orig_conf_taxa)))
for tax in orig_conf_taxa:
    sys.stdout.write("{}\n".format(tax))



conflict = resp_updated.response_dict

updated_conf_taxa = set()
for node in conflict:
    if conflict[node]['status'] == 'conflicts_with':
        witness = conflict[node]['witness_name']
        for tax in witness:
            updated_conf_taxa.add(tax)

sys.stdout.write("\nUpdated tree conflicts with {} taxa in the OpenTree taxonomy:\n".format(len(updated_conf_taxa)))
for tax in updated_conf_taxa:
    sys.stdout.write("{}\n".format(tax))

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

tree_updated.write(path = "{}/conflict_label.tre".format(comparisondir), schema="newick")



## tree_updated_synth = conflict_tree(unpruned_tree2, otu_dict)
## resp_updated_synth = OT.conflict_str(tree_updated_synth.as_string(schema="newick"), 'synth')
## conflict_synth = resp_updated_synth.response_dict
##
##
## for node in tree_updated_synth:
##     if node.taxon:
##         node_id = node.taxon.label.split('_')[1]
##         conf_node = conflict_synth.get(node_id, False)
##         if conf_node:
##             new_label = "{} {}".format(conf_node['status'], conf_node['witness'])
##             node.taxon.label = new_label
##     else:
##         node_id = node.label.split('_')[1]
##         conf_node = conflict_synth.get(node_id, False)
##         if conf_node:
##             new_label = "{} {}".format(conf_node['status'], conf_node['witness'])
##             node.label = new_label
##
