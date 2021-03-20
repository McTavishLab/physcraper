## Analysing the Physcraper results

*Under construction*

The functionalities described in the following sections are still under development.
Even in this raw form, we think they have the potential to be pretty useful, so
we released them.

The first two functionalities are not accessible through the command line yet.
This means that you can access them only through Python for the moment.
The last functionality can be accessed through the command line, but it has not
achieved final form yet.

### Relabeling the trees

For downstream analyses and figure making, it can be handy to swap name on tips of
the updated phylogeny from alternative taxonomies or taxon id numbers.

Do that with the `write_labelled` function, e.g.,

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    pg55.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')

### Rerooting the trees

A correctly rooted phylogeny is needed to compare relationships between two or more phylogenies.
Rooting phylogenies can be tricky. Physcraper places a suggested root based on the taxonomic relationships in OpenTree, using the `root_tree_from_synth` function.

To root a Physcraper tree using either the OpenTree taxonomy, or the OpenTree synthetic tree.
First load the tree object:

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')

Then, to root based on taxonomy, set `base = "ott"`:

    from physcraper import opentree_helpers
    ott_rooted_tree = opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='ott')
    pg55.tre = ott_rooted_tree
    pg55.write_labelled(label="^ot:ottTaxonName", path="ott_root.tre")


And, to root based on phylogenetic relationships in the synthetic tree, set `base = "synth"`:

    from physcraper import opentree_helpers
    synth_rooted_tree = opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='synth')
    pg55.tre = synth_rooted_tree
    pg55.write_labelled(label="^ot:ottTaxonName", path="synth_root.tre")


In this example both trees are the same even though they use the MRCA of different pairs of taxa, because those MRCA's map to the same node on the output tree.

However rooting based on OTT can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.
So whenever possible, the root should be specified by the user, for example by choosing from tips in the "otu_dictionary" file, e.g.,

    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    outgroup = ['otu376436','otu376444']
    mrca = pg55.tre.mrca(taxon_labels=outgroup)
    pg55.tre.reroot_at_node(mrca, update_bipartitions=True)
    pg55.write_labelled(label="^ot:ottTaxonName", path="manual_root.tre")


### Tree comparison with Robinson-Foulds (RF) distance

The `tree_comparison.py` script takes as an argument the output directory of a Physcraper run,
and compares the relationships in the final tree to the relationships in the input tree.

It uses the rooting functions described above to assure the two trees have the same root.

The simplest command line is:

    tree_comparison.py -d docs/examples/pg_55_web/ -o pg_55_comparison

Where `-d` is the directory with the output files from a Physcraper run, and `-o`
is the path and file where the tree comparison results are stored. By default it will root trees based on the OpenTree taxonomy, as explained above.

Alternatively, you can pass in OTT ids of two taxa from the input tree to use to root both trees,
e.g.,

    tree_comparison.py -d docs/examples/pg_55_web/ -og otu376420 otu376439 otu376452 -o pg_55_comparison

In either case the script will print to screen information comparing the two trees including:

* The number of new tips
* The number of new taxa
* Whether the taxa in the tree are included synthesis phylogenies currently in OpenTree
* Which taxa phylogenetic information is not currently incorporated into the synthetic tree
* The RF distance and weighted RF distance between the relationships of tips that are in both trees
* How the estimates of phylogenetic relationships of taxa included in the OpenTree taxonomy from both trees have conflict with the monophyly of the OpenTree synthetic tree.
