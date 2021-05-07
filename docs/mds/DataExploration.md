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

For downstream analyses and figure making, it can be handy to swap labels on tips of
the updated phylogeny from alternative taxonomies or taxon id numbers.

Do that with the `write_labelled` function. For example, to change  e.g.,

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    pg55.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='tests/tmp/pg_55_repeats.tre')

### Rerooting the trees

A correctly rooted phylogeny is needed to compare relationships between two or more phylogenetic hypotheses.
Automatic rooting of phylogenies is not straightforward. Physcraper's `root_tree_from_synth` function places a suggested root based on relationships in OpenTree's synthetic tree or in its taxonomic tree.

To root a Physcraper tree using either the OpenTree taxonomy, or the OpenTree synthetic tree.
First load the tree object:

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')

Then, to root based on the OpenTree taxonomy, set `base = "ott"`:

    from physcraper import opentree_helpers
    ott_rooted_tree = opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='ott')
    pg55.tre = ott_rooted_tree  # set tree rooted based on taxonomy as the tree object
    pg55.write_labelled(label="^ot:ottTaxonName", path="tests/tmp/pg_55_ott_root.tre")


And, to root based on phylogenetic relationships in the OpenTree synthetic tree, set `base = "synth"`:

    from physcraper import opentree_helpers
    synth_rooted_tree = opentree_helpers.root_tree_from_synth(pg55.tre, pg55.otu_dict, base='synth')
    pg55.tre = synth_rooted_tree  # set tree rooted based on synth as the tree object
    pg55.write_labelled(label="^ot:ottTaxonName", path="tests/tmp/pg_55_synth_root.tre")


In this example both trees are the same even though they use the MRCA of different pairs of taxa, because those MRCA's map to the same node on the output tree.

However rooting based on OTT can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.
So whenever possible, the root should be specified by the user, for example by choosing from tips in the "otu_info.csv" file in the `outputs` folder of a Physcraper run, e.g.,

    pg55 = treetaxon.generate_TreeTax_from_run('docs/examples/pg_55_web')
    outgroup = ['otu376436','otu376444']
    mrca = pg55.tre.mrca(taxon_labels=outgroup)
    pg55.tre.reroot_at_node(mrca, update_bipartitions=True)
    pg55.write_labelled(label="^ot:ottTaxonName", path="tests/tmp/pg_55_manual_root.tre")


### Tree comparison with Robinson-Foulds (RF) distance

The `tree_comparison.py` script takes as an argument the output directory of a Physcraper run,
and compares the relationships in the final tree to the relationships in the input tree.

It uses the rooting functions described above to assure the two trees have the same root.
By default it will root trees based on the OpenTree taxonomy.

Usage:
    tree_comparison.py  [-h] [-d PHYSCRAPER_OUTPUT_DIRECTORY] [-t1 ORIGINAL_TREE] [-t2 UPDATED_TREE]
                        [-otu OTU_INFO_FILE] [-og OUTGROUP] [-o OUTPUT_DIRECTORY]

Arguments:
  -h, --help            show this help message and exit
  -d PHYSCRAPER_OUTPUT_DIRECTORY, --dir PHYSCRAPER_OUTPUT_DIRECTORY
                        A folder name containing all output directories from a Physcraper run
  -o OUTPUT_DIRECTORY, --output OUTPUT_DIRECTORY
                        A folder name to write the results of the comparison analysis to
  -t1 ORIGINAL_TREE, --original_tree ORIGINAL_TREE
                        Original tree from `inputs` folder
  -t2 UPDATED_TREE, --updated_tree UPDATED_TREE
                        Updated tree from `outputs` folder
  -otu OTU_INFO_FILE, --include_missing INCLUDE_MISSING
                        File with taxon information in JSON format from the `run` folder

The simplest comparison run from the command line:

    tree_comparison.py  [-h] [-d PHYSCRAPER_OUTPUT_FOLDER] [-o OUTPUT_DIRECTORY]

This compares the original tree in the `inputs` folder and the updated tree in the `outputs` folder
For example:

    tree_comparison.py -d docs/examples/pg_55_web/ -o pg_55_comparison


Alternatively, you can pass in OTT ids of two or more taxa from the input tree to use as outgroups to root both trees.
For example:

    tree_comparison.py -d docs/examples/pg_55_web/ -og otu376420 otu376439 otu376452 -o pg_55_comparison

If the comparison between the two trees is possible, the script will print to screen information comparing the two trees, including:

* The number of new tips
* The number of new taxa
* Whether the taxa in the tree are included synthesis phylogenies currently in OpenTree
* Which taxa phylogenetic information is not currently incorporated into the synthetic tree
* The RF distance and weighted RF distance between the relationships of tips that are in both trees
* How the estimates of phylogenetic relationships of taxa included in the OpenTree taxonomy from both trees have conflict with the monophyly of the OpenTree synthetic tree.
