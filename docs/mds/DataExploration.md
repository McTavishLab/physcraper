## Analysing the Physcraper results

### Rerooting the trees

*In construction*

A correctly rooted phylogeny is needed to compare relationships between two or more phylogenies.
Rooting phylogenies can be tricky. Physcraper places a suggested root based on the taxonomic relationships in OpenTree, using this function.

To root a physcraper tree using either the OpenTree taxonomy, or the OpenTree synth tree, use `root_tree_from_synth`:
e.g.
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    # To root based on taxonomy, set base = "ott"
    ott_rooted_tree = physcraper.opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='ott')
    pg55.tre = ott_rooted_tree
    pg55.write_labelled(label="^ot:ottTaxonName", path="ott_root.tre")


    # To root based on phylogenetic relationships in the synthetic tree, set base = "synth"
    synth_rooted_tree = physcraper.opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='synth')
    pg55.tre = synth_rooted_tree
    pg55.write_labelled(label="^ot:ottTaxonName", path="synth_root.tre")


In this example both trees are the same even though they use the MRCA of differnet pairs of taxa, because those MRCA's map to the same node on the output tree.


However rooting based on ott can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.
So whenever possible, the root should be specified by the user, for example by choosing from tips in the otu_dictionary file.

e.g 
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    outgroup = ['otu376436','otu376444']
    mrca = pg55.tre.mrca(taxon_labels=outgroup)
    pg55.tre.reroot_at_node(mrca, update_bipartitions=True)
    pg55.write_labelled(label="^ot:ottTaxonName", path="manual_root.tre")


### Tree comparison with Robinson-Foulds

*In construction*

The 'tree_comparison.py' script takes as an argument the output directory of a physcraper run, 
and compares the relationships in teh final tree to the relationships in the input tree.

It uses the rooting functions descibed above to assure the two trees are rooted the same.

the simplest format is

    tree_comparison.py -d docs/examples/pg_55_web/ -o pg_55_comparison

Where -d is the directory with the output files from a physcraper run, and -o is where the tree compaison results are stored. By default it will root trees based on the OpenTree taxonomy.

Alternatively, you can pass in tip OTU ids from the input tree to use to root both trees:
e.g.

  tree_comparison.py -d docs/examples/pg_55_web/ -og otu376420 otu376439 otu376452 -o pg_55_comparison

In either case the script will print to screen information comparing the two trees including:
    - The number of new tips
    - The number of new taxa
    - Whether the taxa in the tree are included synthesis phylogenies currently in OpenTree
    - Which taxa phylogenetic information is not currently incorporated into the synthetic tree
    - The RF distance and weighted RF distance between the relationships of tips that are in both trees
    - Which taxa included in the OpenTree taxonomy these estimates conflict with the monophyly of.


### Relabeling the trees

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    pg55.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')
