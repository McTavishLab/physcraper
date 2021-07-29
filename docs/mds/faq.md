Frequently asked questions

## I have to learn to use OpenTree to use Physcraper, is the learning curve worth it?

We think that this decision depends on the goals of the user.

To our knowledge, existing tools that automatize dataset construction for phylogenetics focus on assembling an alignment _de novo_ and/or mining and filtering homolog sequences.
The main goal of Physcraper is to construct upon the knowledge contained in previous expertly-curated and peer-reviwed homology hypothesis (alignments), as well as establish interoperability between biological databases and phylogenetic knowledge. To achieve that, understanding the OpenTree Taxonomy and associated tools is key.

In this sense Physcraper offers the unique advantage to automatically connect a phylogenetic tree to other databases and services, including alignment databases such as TreeBASE, and the conflicting service provided by OpenTree, which permits to rapidly identify regions in the tree that:
- have been enriched in a tree with phylogenetic information,
- are coherent with other phylogenetic estimates, as well as
- conflict with other phylogenetic estimates.
-
Users can update phylogenetic trees using any other existing tool, and if they want to connect them to other biological databases, particularly in a reproducible way, matching their taxon names to the OpenTree Taxonomy allows them to do this programmatically, instead of by hand.

For this goal, having a minimum familiarity with the OpenTree tools is needed.

We realize that this might initially discourage some users, but we believe that the benefits brought by connecting taxonomic data with the OpenTree services will encourage users to familiarize with the OpenTree services, and to adopt the use of Physcraper.


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
