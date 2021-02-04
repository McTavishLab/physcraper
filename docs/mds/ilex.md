---
title: "Physcraper Examples: the Hollies"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ilex}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---





<img src="images/holly-snow.jpg" title="*Ilex* fruits and leaves in the snow" alt="*Ilex* fruits and leaves in the snow" width="75%" />

## I. Finding a tree to update

### With the Open Tree of Life website

Go to the [Open Tree of Life website](https://tree.opentreeoflife.org/opentree/argus/opentree12.3@ott93302) and use the "search for taxon" menu to look up the taxon *Ilex*.

This is how the genus *Ilex* is represented on the Open Tree of Life synthetic tree at the middle of year 2020:

<img src="ilex-tree-opentree-2020-05-31.png" title="Go to the website and look at it [here](https://tree.opentreeoflife.org/opentree/argus/ottol@727571/Ilex)" alt="Go to the website and look at it [here](https://tree.opentreeoflife.org/opentree/argus/ottol@727571/Ilex)" width="100%" />

***

Navigating into the tree, we notice that there might be two studies associated to this portion of the Open Tree synthetic tree.

<img src="ilex-tree-opentree-2020-05-31-2.png" width="100%" />

***

Let's verify that on the [study curator of OToL](https://tree.opentreeoflife.org/curator).

Studies matching the word 'ilex' on the curator database, at the middle of year 2020. Some of these studies are not actually about the hollies, but other taxa that have the species epithet *ilex*.

<img src="ilex-curation-opentree-2020-05-31.png" width="100%" />

### Finding a tree to update using the R package rotl

<!-- ```{r} -->
<!-- ilex <- rotl::tnrs_match_names(names = "ilex") -->
<!-- ``` -->

<!-- ```{r, results='asis', echo=FALSE} -->
<!-- knitr::kable(ilex, caption = "Taxonomic Name Resolution information for the genus *Ilex*") -->
<!-- ``` -->

<!-- Can we use the *OTT id* to get some trees or studies? -->

<!-- ```{r} -->
<!-- rotl::studies_properties() -->
<!-- ``` -->

Explain what a focal clade is.

There is a handy function that will search a taxon among the focal clades reported across trees.


```r
studies <- rotl::studies_find_studies(property="ot:focalCladeOTTTaxonName", value="Ilex")
```


Table: Studies with the genus *Ilex* as focal clade.

|study_ids |n_trees |tree_ids           |candidate |study_year |title                                                                                                                    |study_doi                              |
|:---------|:-------|:------------------|:---------|:----------|:------------------------------------------------------------------------------------------------------------------------|:--------------------------------------|
|ot_1984   |1       |tree1              |          |2020       |                                                                                                                         |http://dx.doi.org/10.1111/jse.12567    |
|pg_2827   |2       |tree6576, tree6577 |tree6577  |2003       |Molecular analyses of the genus Ilex (Aquifoliaceae) in southern South America, evidence from AFLP and ITS sequence data |http://dx.doi.org/10.3732/ajb.92.2.352 |



It seems like the oldest tree, *tree6577* from study *pg_2827*, is in the Open Tree of Life synthetic tree.

Let's get it and plot it here:



```r
original_tree <- rotl::get_study_tree(study_id = "pg_2827", tree_id = "tree6577")
#> Warning in build_raw_phylo(ncl, missing_edge_length): missing edge lengths are
#> not allowed in phylo class. All removed.
ape::plot.phylo(ape::ladderize(original_tree), type = "phylogram", cex = 0.3, label.offset = 1, edge.width = 0.5)
```

<img src="../mds/ilex-images/original-tree-1.png" title="plot of chunk original-tree" alt="plot of chunk original-tree" width="100%" style="background-color: #DCDCDC; padding:10px; display: inline-block;" />

***

Now, let's look at some properties of the tree:


```r
ape::Ntip(original_tree)  # gets the number of tips
#> [1] 48
ape::is.rooted(original_tree)  # check that it is rooted
#> [1] TRUE
ape::is.binary(original_tree)  # check that it is fully resolved
#> [1] FALSE
datelife::phylo_has_brlen(original_tree)  # checks that it has branch lengths
#> [1] FALSE
```
The tree has 48 tips, is rooted, has no branch lengths and is not fully resolved, as you probably already noted. Also, labels correspond to the labels reported on the original study [here](http://dx.doi.org/10.3732/ajb.92.2.352). Other labels are available to use as tip labels. For example, you can plot the tree using the unified taxonomic names, or the taxonomic ids.

### Finding a tree in the Open Tree of Life phylesystem with at least one tip label belonging to the group of interest.

```
find_trees.py --taxon_name  "Ilex"
```

## II. Getting the underlying alignment

### TreeBASE

#### Using physcraper and the arguments `-tb` and `-no_est`

```
physcraper_run.py -s pg_2827 -t tree6577 -tb -no_est -o data/pg_2827_tree6577
```

#### Downloading the alignment directly from a repository

The alignments are here <https://treebase.org/treebase-web/search/study/matrices.html?id=1091>

On a mac you can do:
```
wget "http://purl.org/phylo/treebase/phylows/matrix/TB2:M2478?format=nexus" -o data-raw/alignments/T1281-M2478.nex
```


### Other repos


## III. Starting a physcraper a run


```
physcraper_run.py -s pg_2827 -t tree6577 -o data/pg_2827_tree6577
```

```
physcraper_run.py -s pg_2827 -t tree6577 -a data-raw/alignments/T1281-M2478.nex -as nexus -o data/ilex-remote
```

### Using a local BLAST database

```
physcraper_run.py -s pg_2827 -t tree6577 -a data-raw/alignments/T1281-M2478.nex -as nexus -db local_blast_db -o data/ilex-local
```

## IV. Reading physcraper results

### The physcraper tag

T1281-M2478

### Input files

Physcraper rewrites input files for a couple reasons: reproducibility, taxon name matching, and taxon reconciliation.
It writes the config file down if none was provided

### Run files

Files in here are also automatically written down by physcraper.

blast runs, alignments, raxml trees, bootstrap

The trees are reconstructed using RAxML, with tip labels corresponding to local ids (e.g., otu42009, otuPS1) and not taxon names (e.g., *Helwingia japonica*), nor taxonomic ids (e.g., ott or ncbi). Branch lengths are proportional to relative substitution rates.
The RAxML tree with taxon names as tip labels is saved on the `outputs_tag` folder.



```r
updated_tree_otus <- ape::read.tree(file = "../data/ilex-local/run_T1281-M2478/RAxML_bestTree.2020-06-29")
ape::plot.phylo(ape::ladderize(updated_tree_otus), type = "fan", cex = 0.25, label.offset = 0.01, edge.width = 0.5)
ape::add.scale.bar(cex = 0.3, font = 1, col = "black")
mtext("Updated tree - otu tags as labels", side = 3, line = -1)
```

<img src="../mds/ilex-images/updated-tree-otus-1.png" title="plot of chunk updated-tree-otus" alt="plot of chunk updated-tree-otus" width="100%" style="background-color: #6B8E23; padding:10px; display: inline-block;" />

### Output files

A number of files is automatically written down by physcraper.

A nexson tree with all types of tip labels is saved in here.
From this tree, a tree with any kind of label can be produced.
By default, the updated tree with taxon names as tip labels is saved
in the `output_tag` folder as `updated_taxonname.tre`.



```r
updated_tree_taxonname <- ape::read.tree(
  file = "../data/ilex-local/outputs_T1281-M2478/updated_taxonname.tre")
# ape::plot.phylo(ape::ladderize(updated_tree_taxonname), cex = 0.35)
```


```r
ape::plot.phylo(ape::ladderize(updated_tree_taxonname), type = "fan", cex = 0.25, label.offset = 0.01, edge.width = 0.5)
ape::add.scale.bar(cex = 0.3, font = 1, col = "black")
mtext("Updated tree - Taxon names as labels  ", side = 3, line = -1)
```

<img align = "center" width="100%" src="https://raw.githubusercontent.com/McTavishLab/physcraperex/master/vignettes/mds/ilex-images/updated-tree-taxonname-1.png" title="plot of chunk updated-tree-taxonname" alt="plot of chunk updated-tree-taxonname" style="background-color: #9ecff7; padding:10px; display: inline-block;">



## IV. Analyzing the physcraper results

First, compare the original tree with the updated tree.

We can prune the updated tree, so it is a more straight forward comparison:


```r
original_tree_taxonname <- ape::read.tree(file = "../data/ilex-local/inputs_T1281-M2478/taxonname.tre")
tips2keep <- match(original_tree_taxonname$tip.label,updated_tree_taxonname$tip.label)

pruned_updated_tree_taxonname <- ape::drop.tip(updated_tree_taxonname, updated_tree_taxonname$tip.label[-tips2keep])
```


```r
original_tree_taxonname$edge.length <- NULL
original_tree_taxonname2 <- original_tree_taxonname
original_tree_taxonname2$tip.label <- gsub("_otu.*", "", original_tree_taxonname$tip.label)

pruned_updated_tree_taxonname2 <- pruned_updated_tree_taxonname
pruned_updated_tree_taxonname2$tip.label <- gsub("_otu.*", "", pruned_updated_tree_taxonname$tip.label)

cotree1 <- phytools::cophylo(ape::ladderize(original_tree_taxonname2), ape::ladderize(pruned_updated_tree_taxonname2), rotate.multi =TRUE)
```

Rotating nodes to optimize matching...
Done.

```r

phytools::plot.cophylo(
  x = cotree1,
  fsize = 2,
  lwd = 0.4,
  mar=c(.1,.1,8,.3),
  ylim=c(-.1,1),
  scale.bar=c(0, 0.5),
  link.lwd=5,
  link.lty="solid",
  link.col=phytools::make.transparent("#8B008B",0.5))
mtext(expression(bold("Original tree")), cex = 4.5, side= 3, line = 0, adj = 0.2)
mtext(expression(bold("Pruned updated tree")), cex = 4.5, side= 3, line = 0, adj = 0.9)
```

<img src="../mds/ilex-images/cotree-plot1-1.png" title="plot of chunk cotree-plot1" alt="plot of chunk cotree-plot1" width="90%" style="background-color: #8B008B; padding:10px; display: inline-block;" />


But, comparing vs the whole updated tree is more interesting:


```r
updated_tree_taxonname2 <- updated_tree_taxonname
updated_tree_taxonname2$tip.label <- gsub("_otu.*", "", updated_tree_taxonname$tip.label)
# updated_tree_taxonname2$edge.length <- updated_tree_taxonname$edge.length*10 # for viz purposes! did not work
# updated_tree_taxonname2$edge.length <- NULL # does not work either
cotree2 <-  phytools::cophylo(
  ape::ladderize(original_tree_taxonname2),
  ape::ladderize(updated_tree_taxonname2)) #,
```

Rotating nodes to optimize matching...
Done.

```r
  #rotate.multi =TRUE) # it does not work with these trees for some reason.
phytools::plot.cophylo(
  x = cotree2,
  fsize = 1,
  lwd = 0.4,
  mar=c(.1,.1,8,.5),
  ylim=c(-.1,1),
  # xlim=c(-1,1),
  scale.bar=c(0, 1),
  # link.type="curved",
  link.lwd=2,
  link.lty="solid",
  link.col=phytools::make.transparent("#8B008B",0.5))
mtext(expression(bold("Original tree")), cex = 4, side = 3, line = 0, adj =0.2)
mtext(expression(bold("Updated tree")), cex = 4, side = 3, line = 0, adj = 0.8)
```

<img src="../mds/ilex-images/cotree-plot2-1.png" title="plot of chunk cotree-plot2" alt="plot of chunk cotree-plot2" width="100%" style="background-color: #8B008B; padding:10px; display: inline-block;" style="display: block; margin: auto;" />

### Now, compare the updated tree with the newer tree


```r
newer_tree <- rotl::get_study_tree(study_id = "ot_1984", tree_id = "tree1")
ape::plot.phylo(ape::ladderize(newer_tree), type = "fan", cex = 0.3, label.offset = 0.5, edge.width = 0.5)
mtext("Newer tree - matched labels", side = 3, cex = 0.6)
```

<img src="../mds/ilex-images/newer-tree-1.png" title="plot of chunk newer-tree" alt="plot of chunk newer-tree" width="100%" style="background-color: #FF1493; padding:10px; display: inline-block;" />


```r

cotree3 <-  phytools::cophylo(
  ape::ladderize(newer_tree),
  ape::ladderize(updated_tree_taxonname2)) #,
```

Rotating nodes to optimize matching...
Done.

```r
  #rotate.multi =TRUE) # it does not work with these trees for some reason.
phytools::plot.cophylo(
  x = cotree3,
  fsize = 1,
  lwd = 0.4,
  mar=c(.1,.1,8,.5),
  ylim=c(-.1,1),
  # xlim=c(-1,1),
  scale.bar=c(0, 1),
  # link.type="curved",
  link.lwd=2,
  link.lty="solid",
  link.col=phytools::make.transparent("#FF1493",0.5))
mtext(expression(bold("Newer tree")), cex = 4, side = 3, line = 0, adj =0.2)
mtext(expression(bold("Updated old tree")), cex = 4, side = 3, line = 0, adj = 0.8)
```

<img src="../mds/ilex-images/cotree-plot3-1.png" title="plot of chunk cotree-plot3" alt="plot of chunk cotree-plot3" width="100%" style="background-color: #FF1493; padding:10px; display: inline-block;" style="display: block; margin: auto;" />

## V. Citation

## VI. Acknowledgments

University of California, Merced cluster, MERCED (Multi-Environment Research Computer for Exploration and Discovery) supported by the National Science Foundation (Grant No. ACI-1429783).

Holly Image by <a href="https://pixabay.com/es/users/WolfBlur-2503887/?utm_source=link-attribution&amp;utm_medium=referral&amp;utm_campaign=image&amp;utm_content=3012084">Wolfgang Claussen</a> at <a href="https://pixabay.com/es/?utm_source=link-attribution&amp;utm_medium=referral&amp;utm_campaign=image&amp;utm_content=3012084">Pixabay</a>
