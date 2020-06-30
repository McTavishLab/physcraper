### Find a study with your taxon of interest


You can upload your own tree to OpenTree to update it.
See [Submitting-phylogenies-to-Open-Tree-of-Life](https://github.com/OpenTreeOfLife/opentree/wiki/Submitting-phylogenies-to-Open-Tree-of-Life)


For this example we'll use find a tree that is already in the database.
Search on OpenTree of life using your taxon of interest, e.g. 'Malvaceae'

    $ find_trees.py --taxon_name "Malvaceae"


Lets use Wilkie et al (2006). https://tree.opentreeoflife.org/curator/study/view/pg_55

While this study was focussed on the family "Sterculiacea",
phylogenetic inference have suggested that this family is not monophyletic.

https://tree.opentreeoflife.org/opentree/argus/ottol@996482

In order to further assess, lets update the tree!


### Run the auto update


The script `physcraper_run.py` wraps together linking the tree and alignment, blasting, aligning sequences, and inferring an updated tree.
Detailed explanation of the inputs needed can be found at [PhyscraperRun](./PhyscraperRun.md).

The blast search part of updating trees takes a long time (for example, this analysis took around 12 hours!).


    $ physcraper_run.py -s pg_55 -t tree5864 -tb -r -o pg_55


We have put example outputs from this command in `docs/examples/pg_55`, so that you can explore the outputs without waiting for the searches to complete.

### Output stucture

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
   
   -- CSV file with information about sequences




### Compare your new tree to existing relationships

A correctly rooted phylogeny is needed to compare taxonomic groups.
Rooting phylogenies can be tricky. While physcraper places a suggested root based on the taxonomic relationships in OpenTree, it often is unreliable.

There is a simple tree comparison script at `tree_comparison.py`

Detailed explanation of that script, amre more ways to explore the data are described in [DataExploration](./DataExploration.md)


    tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o pg_55_comparison
