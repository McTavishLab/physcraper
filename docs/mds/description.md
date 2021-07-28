## The Physcraper framework

While genome scale data is increasing rapidly, there are still large quantities
of gene-sequence data being uploaded to the US National Center on Biotechnology
Information (NCBI) database [GenBank](https://www.ncbi.nlm.nih.gov/genbank/statistics/).
These data are often appropriate for looking at phylogenetic relationships, and
have the advantage of being homologous to genetic sequences used to construct existing
trees.

If you have access to a single gene DNA alignment and a tree, Physcraper automates
adding new lineage samples into your tree by using [Open Tree of Life](#the-open-tree-of-life) tools coupled to the [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) algorithm to search for loci in [GenBank](https://www.ncbi.nlm.nih.gov/genbank/statistics/) that are likely to be homologous to sequences in the initial DNA alignment.

<br/>

![](../img/schematic-final.svg)

Figure 1 from [Sanchez-Reyes et al. 2021](https://doi.org/10.1186/s12859-021-04274-6):
The Physcraper framework consists of 4 general steps. The methodology is extensively described in the [Implementation](https://physcraper.readthedocs.io/en/latest/implementation.html) section of this documentation.

<br/>

By using a starting alignment and tree, Physcraper takes advantage of loci that
previous researchers have assessed and deemed appropriate for the phylogenetic scope.
The sequences added in the search are limited either to a user specified taxon or
monophyletic group, or within the taxonomic scope of the ingroup of the starting tree.

These automated trees can provide a quick inference or potential relationships,
of problems in the taxonomic assignments of sequences, and flag areas of potential systematic interest.

<br/>


## The Open Tree of Life

The Open Tree of Life ([OpenTree](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302)) is a project that unites phylogenetic inferences and taxonomy
to provide a synthetic estimate of species relationships across the entire tree of life.

<br/>

![](../img/synthtreeleg.svg)

OpenTree synthetic tree. Figure 1 from [Hinchliff et al. 2015](https://www.pnas.org/content/112/41/12764.short).
For more information on the OpenTree project go to https://opentreeoflife.github.io

<br/>

OpenTree aims to construct a comprehensive, dynamic and digitally-available tree
of life by synthesizing published phylogenetic trees along with taxonomic data.
Currently the tree comprises 2.3 million tips.
However, only around 90,000 of those taxa are represented by phylogenetic estimates -
the rest are placed in the tree based on their taxonomic names.

To achieve this, the OpenTree Taxonomy (OTT) constructs a reference taxonomy through
an algorithmic combination of several source taxonomies, such as:
- [Hibbet et al. 2007](https://doi.org/10.1016/j.mycres.2007.03.004),
- [SILVA](http://www.arb-silva.de/),
- the [Index Fungorum](http://www.indexfungorum.org/),
- [Schäferhoff et al. 2010](https://doi.org/10.1186/1471-2148-10-352),
- the [World Register of Marine Species](WoRMS; http://www.marinespecies.org/aphia.php)
- the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/books/NBK21100/),
- the Global Biodiversity Information facility [(GBIF) backbone Taxonomy](https://www.gbif.org/), and
- the [Interim Register of Marine and Nonmarine Genera (IRMNG)](https://irmng.org/).


<br/>
