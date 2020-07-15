### Compare your new tree to existing relationships

A correctly rooted phylogeny is needed to compare taxonomic groups.
Rooting phylogenies can be tricky. While physcraper places a suggested root based on the taxonomic relationships in OpenTree,
this root can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.

There is a simple tree comparison script, `tree_comparison.py`

Detailed explanation of that script, and more ways to explore the data are described in [data exploration](https://physcraper.readthedocs.io/en/latest/data_exploration.html)


    tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o pg_55_comparison


The taxonomic name mapping makes comparisons across trees straightforward.  
