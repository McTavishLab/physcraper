Frequently asked questions

## How does Physcraper handle polytomies of starting trees?

The Physcraper starting tree is a phylogeny whose tip labels must have been standardized to the OpenTree Taxonomy (as described in the Introduction section:
[Mapping names to taxa](https://physcraper.readthedocs.io/en/latest/quick-start.html#updating-your-own-tree-and-alignment)).
Original tip labels of the starting tree must be identical to taxon labels on the starting alignment.
However, not all taxon labels in the alignment have to be present in the tree and
visceversa.

Physcraper makes use of the starting tree in four main ways:
1. to delimit a taxon for the GenBank search (a search taxon),
2. to be used as starting tree for the phylogenetic reconstruction software of choice,
3. to standardize the taxon names from the starting alignment, and
4. to compare the updated phylogenetic relationships with the original ones.


Physcraper does not really "handle" polytomies. The goal of the software is to use the
existing phylogenetic information that has been generated, reviewed, published and curated by experts in the field.

If a starting tree contains polytomies, these can only affect the outcome of the analysis if the starting tree is used for the case (1) delimiting a taxon for the GenBank search.
To delimit the search taxon from the starting tree, a known outgroup is necessary.
The outgroup can be user defined. If the outgroup is not defined by the user, Physcraper will attempt to root the starting tree following the OpenTree Taxonomy.
If succesful, it will take the tip labels from the earliest diverging branch with the least number of tips. These will be used as outgroup. However, if the starting tree has polytomies
around the early diverging branches,
the automatic rooting is problematic and can have multiple solutions.



## How does Physcraper use the starting alignment?

Physcraper uses the input DNA alignment (single or multiple marker) to mine the GenBank database with the goal of increasing the
lineage sampling of the alignment within a given biological group.
