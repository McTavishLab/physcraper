import dendropy
from dendropy import Tree
from dendropy.calculate import treecompare

place_tree_path = "kwyj_copy/rpb2/RAxML_entropy.PLACE"
full_tree_path = "kwyj_copy/rpb2/random_resolve2016-02-04.tre"


original_tree_path = "kwyj_copy/rpb2/original_pruned.tre"
original_recalc_path = "kwyj_copy/rpb2/RAxML_bestTree.orig"

full_tree = Tree.get(path=full_tree_path,
                     schema="newick",
                     preserve_underscores=True)

original_tree = Tree.get(path=original_tree_path,
                         schema="newick",
                         preserve_underscores=True,
                         taxon_namespace=full_tree.taxon_namespace)

orig_taxa = []
for node in original_tree.leaf_nodes():
   orig_taxa.append(node.taxon)

full_pruned = Tree.get(path=full_tree_path,
                       schema="newick",
                       preserve_underscores=True,
                       taxon_namespace=full_tree.taxon_namespace)

full_pruned.retain_taxa(orig_taxa)
full_pruned.write(path="kwyj_copy/rpb2/full_pruned.tre",
                  schema="newick",
                  unquoted_underscores=True)


original_tree.encode_bipartitions()
full_pruned.encode_bipartitions()
print("now for orig vs pruned")
print(treecompare.weighted_robinson_foulds_distance(original_tree, full_pruned))
print(treecompare.false_positives_and_negatives(original_tree, full_pruned))

print("now for orig recalc vs pruned")

original_recalc = Tree.get(path=original_recalc_path,
                         schema="newick",
                         preserve_underscores=True,
                         taxon_namespace=full_tree.taxon_namespace)

print(treecompare.weighted_robinson_foulds_distance(original_recalc, full_pruned))
print(treecompare.false_positives_and_negatives(original_recalc, full_pruned))

print("now for orig vs orig recalc")

print(treecompare.weighted_robinson_foulds_distance(original_tree, original_recalc))
print(treecompare.false_positives_and_negatives(original_tree, original_recalc))

#full

#self.tre.prune_taxa_with_labels(prune)


#raxmlHPC -m GTRCAT -s original_pruned.fas -t original_pruned.tre -p 1 -n orig