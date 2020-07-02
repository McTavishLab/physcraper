
### Updating gene trees

If you have access to a single gene alignment, and a tree, you can automate adding homologous data into your tree by searching GenBank.

While genome scale data is increasing rapidly - there are still large quantities of gene-sequence data being uploaded to NCBI GenBank.

<img src="img/seq_data.png" alt="drawing" width="400"/>  


These data are often appropriate for looking at phylogenetic relationships.

Using Physcraper we can use Blast to search for loci that are likely to be homologous to sequences in an existing alignment.

By using a starting tree and alignment, Physcraper, takes advantage of loci that previous researchers have assessed and deemed appropriate for the phylogentic scope.
The sequences added in the search are limited to a user specified taxon or monophyletic group, or within the taxonomic scope of the in-group of the starting tree.

These automated tree can provide a quick inference or potential relationships, of problems in the taxonomic assignments of sequences, and flag areas of potential systematic interest.


## The Open Tree of Life

The Open Tree of Life (https://opentreeoflife.github.io/) is a project that unites phylogenetic inferences and taxonomy to provide a synthetic estimate of species relationships across the entire tree of life.  
![](img/otol_logo.png)  


Open Tree of Life aims to construct a comprehensive, dynamic and digitally-available tree of life by synthesizing published phylogenetic trees along with taxonomic data.  
Currently the tree comprises 2.3 million tips. 
However, only around 90,000 of those taxa are represented by phylogenetic estimates - the rest are placed in the tree based on their taxonomic names.

https://opentreeoflife.github.io/browse/



## Updating a tree from OpenTree of Life

The Open Tree of Life data store, [Phylesystem](https://academic.oup.com/bioinformatics/article/31/17/2794/183373), contains more than 4,500 phylogenetic trees from published studies.  
The tips in these trees are mapped a unified taxonomy, which makes these data searchable in a phylogenetically explicit way.
This is a great place to start of finding existing estimates of phylogenetic relationships, 
and assessing regions of the tree of life which are lacking available phylogenetic estimates.
There is a lot of sequence data that has been generated, but has never been incorporated into any phylogenetic estimates.


### Find a study with your taxon of interest


For this example we'll use find a tree that is already in the OpenTree of Life database.
Search on OpenTree of life using your taxon of interest, e.g. 'Malvaceae'

    $ find_trees.py --taxon_name "Malvaceae"

This prints a bunch of studies out to the screen. We will need an alignment to update (which OpenTree doesn't store), so lets just look at trees that have data stored in tree base.

    $ find_trees.py --taxon_name "Malvaceae" --treebase

There are a bunch of options!

Lets update Wilkie et al (2006). 
You can view the study on the OpenTree database: [Wilkie2006](https://tree.opentreeoflife.org/curator/study/view/pg_55)

While this study was focussed on the family "Sterculiacea",
phylogenetic inference have suggested that this family is not [monophyletic]((https://tree.opentreeoflife.org/opentree/argus/ottol@996482))

Lets take a look at how recent data affect our inferences of relationships, and if there is sequence data for taxa that don;t have any phylogenetic information available in the tree.

### Run the auto update


The script `physcraper_run.py` wraps together linking the tree and alignment, blasting, aligning sequences, and inferring an updated tree.
Detailed explanation of the inputs needed can be found at [PhyscraperRun](./PhyscraperRun.md).

The blast search part of updating trees takes a long time (for example, this analysis took around 12 hours!).


    $ physcraper_run.py -s pg_55 -t tree5864 -tb -r -o pg_55


We have put example outputs from this command in `docs/examples/pg_55`, so that you can explore the outputs without waiting for the searches to complete.

### Output files

The analysis folder has several sub directories.
each folder is labeled with a 'tag', which by default is the alignment name, but can be set in the `physcraper_run.py` arguments.

The structure consists of:

-  inputs
    -- original tree and alignment

    -- the mapping of the labels to taxa saved as `otu_info.csv`

-  blast_run
    -- blast results for each tip in both the tree and the alignment

-  run
   -- This is where intermediate processing files, and the json formatted otu information are stored

- outputs
   -- final tree and alignment
   
   -- CSV file with information about each sequence



### Compare your new tree to existing relationships

A correctly rooted phylogeny is needed to compare taxonomic groups.
Rooting phylogenies can be tricky. While physcraper places a suggested root based on the taxonomic relationships in OpenTree, 
this root can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.

There is a simple tree comparison script, `tree_comparison.py`

Detailed explanation of that script, and more ways to explore the data are described in [DataExploration](./DataExploration.md)


    tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o pg_55_comparison


## Using your own tree and alignment

You can upload your own tree to OpenTree to update it, and that way it will be included in the synthetic tree!
See [Submitting-phylogenies-to-Open-Tree-of-Life](https://github.com/OpenTreeOfLife/opentree/wiki/Submitting-phylogenies-to-Open-Tree-of-Life)

If you aren't ready to share your tree publicly, you can update it without posting to OpenTree.

You need an alignment (single locus) and a tree. The taxon labels in these two files should be the same.  

You also need a file linking the labels in your tree and alignment to broader taxonomy. This can be easily generated vis OPenTrees Bulk Taxonomic Name Resolution Service. [Bulk TNRS](https://tree.opentreeoflife.org/curator/tnrs/)

Example names:

    physcraper_run.py -tf tests/data/tiny_test_example/test.tre -tfs newick -a tests/data/tiny_test_example/test.fas --taxon_info tests/data/tiny_test_example/main.json -as fasta -o owndata
