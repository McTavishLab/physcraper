## Analysing the Physcraper results

*Under construction*

The functionalities described in the following sections are still under development.
Even in this beta form, we think they have the potential to be useful, so we decided to document them here.

The first two functionalities are not accessible through the command line yet.
This means that you can access them only through Python for the moment.
The last functionality can be accessed through the command line, but it has not
been fully tested yet. Use with caution.

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

Usage:

    tree_comparison.py  [-h] [-d DIRECTORY_NAME] [-t1 FILE_NAME] [-t2 FILE_NAME] [-otu FILE_NAME] [-og OUTGROUP] [-o DIRECTORY_NAME]

Arguments:

<blockquote>
<div><dl class="option-list">
<dt><kbd><span class="option">-h </span>, <span class="option">--help </span></kbd></dt>
<dd><p>Show the help message and exit.</p>
</dd>
<dt><kbd><span class="option">-d <var>DIRECTORY_NAME</var></span>, <span class="option">--dir <var>DIRECTORY_NAME</var></span></kbd></dt>
<dd><p>A name (and path) to a Physcraper output directory.</p>
</dd>
<dt><kbd><span class="option">-o <var>DIRECTORY_NAME</var></span>, <span class="option">--output <var>DIRECTORY_NAME</var></span></kbd></dt>
<dd><p>Name (and path) for a directory to write the results of the comparison analysis to. If it exists, it will be overwritten.</p>
</dd>
<dt><kbd><span class="option">-otu <var>FILE_NAME</var></span>, <span class="option">--otu_info <var>FILE_NAME</var></span></kbd></dt>
<dd><p>Name (and path) of JSON file containing taxon information. Stored in the Physcraper output directory <code class="docutils literal notranslate"><span class="pre">run</span></code> folder.</p>
</dd>
<dt><kbd><span class="option">-t1 <var>FILE_NAME</var></span>, <span class="option">--tree1 <var>FILE_NAME</var></span></kbd></dt>
<dd><p>File name (and path) of original tree from a Physcraper output directory <code class="docutils literal notranslate"><span class="pre">inputs</span></code> folder. Alternatively, any file containing a tree to compare to -t2.</p>
</dd>
<dt><kbd><span class="option">-t2 <var>FILE_NAME</var></span>, <span class="option">--tree2 <var>FILE_NAME</var></span></kbd></dt>
<dd><p>File name (and path) of an updated tree from a Physcraper output directory <code class="docutils literal notranslate"><span class="pre">outputs</span></code> folder. Alternatively, any file containing a tree to compare to -t1.</p>
</dd>
</dl>
</div></blockquote>


This is the simplest command line for comparison of two trees:

    tree_comparison.py  [-h] [-d DIRECTORY_NAME] [-o DIRECTORY_NAME]

For example:

    tree_comparison.py -d docs/examples/pg_55_web/ -o pg_55_comparison

It compares the original tree from the `inputs` folder and the updated tree from the `outputs` folder.
It uses the rooting functions described above to ensure the two trees have the same root.
By default, it will root trees based on the [OpenTree Taxonomy](https://tree.opentreeoflife.org/about/taxonomy-version/ott3.2).

Alternatively, you can pass in OpenTree taxonomic ids (OTT ids) of two or more taxa from the input tree to use as outgroups to root both trees.

For example:

    tree_comparison.py -d docs/examples/pg_55_web/ -og otu376420 otu376439 otu376452 -o pg_55_comparison

If the comparison between the two trees is possible (outgroup-wise), the script will print the results to screen, including:

* The number of new tips
* The number of new taxa
* Whether the taxa in the tree are included synthesis phylogenies currently in OpenTree
* Which taxa phylogenetic information is not currently incorporated into the synthetic tree
* The RF distance and weighted RF distance between the relationships of tips that are in both trees
* How the estimates of phylogenetic relationships of taxa included in the OpenTree taxonomy from both trees have conflict with the monophyly of the OpenTree synthetic tree.
