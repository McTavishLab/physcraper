import sys
import os
import json
import argparse
import dendropy
import copy
import physcraper
from opentree import OT
from dendropy.calculate import treecompare


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

treestr = write_conflict_tree(unpruned_tree2, otu_dict)


resp = OT.conflict_str(treestr, 'ott')

print(resp.response_dict)

#curl -X POST 'https://api.opentreeoflife.org/v3/conflict/conflict-status' 
#-d '{"tree1newick":"((('_nd3_ott318436':0.0039987231091762904,('_nd5_ott318436':0.002661502175610165,(('_nd8_ott114544':0.0013340236695662932,'_nd9_ott114544':0.0026693191584754605,'_nd10_ott114541':0.002671840407990468,'_nd11_ott114544':0.0026706831311005975,('_nd13_ott114544':0.0013300238194868471,'_nd14_ott688671':1.00000050002909e-06)'_nd12_':2.00000100005818e-06,'_nd15_ott114544':0.0013309288367575759)'_nd7_':0.004002067908182651,(('_nd18_ott318436':0.0013373360809382606,'_nd19_ott318436':0.0015450873265851193)'_nd17_':0.004000986236500289,((('_nd23_ott640718':2.00000100005818e-06,'_nd24_ott640718':0.00267032000532656,('_nd26_ott640718':0.001332519568352937,'_nd27_ott640718':0.0013325197242921294)'_nd25_':0.0013346423342407972)'_nd22_':0.0013296539242075864,'_nd28_ott640718':1.00000050002909e-06)'_nd21_':0.002661304633280725,('_nd30_ott688688':1.00000050002909e-06,'_nd31_ott688688':1.00000050002909e-06)'_nd29_':0.009389376977422134)'_nd20_':1.00000050002909e-06)'_nd16_':0.0013253295076721423)'_nd6_':0.0013402480338439555)'_nd4_':1.00000050002909e-06,'_nd32_ott318436':1.00000050002909e-06)'_nd2_':0.0;'","tree2":"synth"}'
