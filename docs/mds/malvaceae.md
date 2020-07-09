## Visualizing the Physcraper results

There are several tree visualization tools available. For reproducibility purposes,
we will use some handy functions from the R language to see our results.

```r
updated_tree_otus <- ape::read.tree(file = "docs/examples/pg_55/run_pg_55tree5864_ndhf/RAxML_bestTree.2020-06-18")
ape::plot.phylo(ape::ladderize(updated_tree_otus), type = "phylogram", cex = 0.25, label.offset = 0.001, edge.width = 0.5)
ape::add.scale.bar(cex = 0.3, font = 1, col = "black")
mtext("Updated tree - otu tags as labels", side = 3)
```

<img src="https://raw.githubusercontent.com/McTavishLab/physcraperex/master/vignettes/mds/malvaceae-images/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="90%" style="background-color: #6B8E23; padding:10px; display: inline-block;" />


```r
updated_tree_taxonname <- ape::read.tree(
  file = "docs/examples/pg_55/outputs_pg_55tree5864_ndhf/updated_taxonname.tre")
```
```r
ape::plot.phylo(ape::ladderize(updated_tree_taxonname), type = "phylogram", cex = 0.2, label.offset = 0.001, edge.width = 0.5)
ape::add.scale.bar(cex = 0.3, font = 1, col = "black")
mtext("Updated tree - Taxon names as labels", side = 3)
```

<img src="https://raw.githubusercontent.com/McTavishLab/physcraperex/master/vignettes/mds/malvaceae-images/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="90%" style="background-color: #9ecff7; padding:10px; display: inline-block;" />


Compare the original tree with the pruned updated tree

```r
original_tree_otus <- ape::read.tree(file = "docs/examples/pg_55/inputs_pg_55tree5864_ndhf/physcraper_pg_55tree5864_ndhf.tre")
updated_tree_otus_pruned <- ape::read.tree(
  file = "../data/pg_55/pruned_updated.tre"
)
```

Now plot them face to face.

We can prune the updated tree, so it is a straight forward comparison:

```r
cotree <-  phytools::cophylo(original_tree_otus, updated_tree_otus_pruned, rotate.multi =TRUE)
```

Rotating nodes to optimize matching...
Done.

```r
phytools::plot.cophylo(x = cotree, fsize = 0.5, lwd = 0.5, mar=c(.1,.2,2,.3), ylim=c(-.1,1), scale.bar=c(0, .05))
title("Original tr.   Pruned updated tr.", cex = 0.1)
```

<img src="https://raw.githubusercontent.com/McTavishLab/physcraperex/master/vignettes/mds/malvaceae-images/cotree-plot1-1.png" title="plot of chunk cotree-plot1" alt="plot of chunk cotree-plot1" width="90%" style="background-color: #8B008B; padding:10px; display: inline-block;" />

But it is more interesting to plot it with all the new tips, so we see exactly where the new things are:

```r
original_tree_taxonname <- ape::read.tree(file = "docs/examples/pg_55/inputs_pg_55tree5864_ndhf/taxonname.tre")
cotree2 <-  phytools::cophylo(
  original_tree_taxonname,
  updated_tree_taxonname,
  rotate.multi =TRUE)
```

Rotating nodes to optimize matching...
Done.

```r
phytools::plot.cophylo(
  x = cotree2,
  fsize = 0.3,
  lwd = 0.4,
  mar=c(.1,.1,2,.5),
  ylim=c(-.1,1),
  scale.bar=c(0, .05),
  # link.type="curved",
  link.lwd=3,
  link.lty="solid",
  link.col=phytools::make.transparent("#8B008B",0.5))
title("Original tree      Updated tree", cex = 0.1)
```

<img src="https://raw.githubusercontent.com/McTavishLab/physcraperex/master/vignettes/mds/malvaceae-images/cotree-plot2-1.png" title="plot of chunk cotree-plot2" alt="plot of chunk cotree-plot2" width="90%" style="background-color: #8B008B; padding:10px; display: inline-block;" />

We can also plot the updated tree against the synthetic subtree of Malvaceae, to visualize how it updates our current knowledge of the phylogeentic relationships within the family.

However, we are having some trouble with matching the tips right now! A fix will come soon:


```r
tolsubtree <- rotl::tol_subtree(ott_id = 279960)
ape::Ntip(tolsubtree)
#> [1] 5898
grep("Pterygota_alata", tolsubtree$tip.label)
#> [1] 5714
updated_tree_taxonname$tip.label
#>  [1] "Fremontodendron_californicum_otuPS13"
#>  [2] "Quararibea_costaricensis_otuPS38"
#>  [3] "Matisia_cordata_otuPS39"
#>  [4] "Hibiscus_bojerianus_otuPS45"
#>  [5] "Macrostelia_laurina_otuPS29"
#>  [6] "Talipariti_tiliaceum_var._tiliaceum_otuPS48"
#>  [7] "Talipariti_hamabo_otuPS47"
#>  [8] "Papuodendron_lepidotum_otuPS46"
#>  [9] "Cephalohibiscus_peekelii_otuPS34"
#> [10] "Kokia_kauaiensis_otuPS32"
#> [11] "Kokia_drynarioides_otuPS31"
#> [12] "Kokia_cookei_otuPS30"
#> [13] "Ochroma_pyramidale_otuPS40"
#> [14] "Ochroma_pyramidale_otu376430"
#> [15] "Catostemma_fragrans_otuPS37"
#> [16] "Scleronema_micranthum_otuPS44"
#> [17] "Cavanillesia_platanifolia_otuPS50"
#> [18] "Spirotheca_rosea_otuPS42"
#> [19] "Spirotheca_rosea_otu376452"
#> [20] "Bombax_buonopozense_otuPS49"
#> [21] "Bombax_buonopozense_otu376420"
#> [22] "Ceiba_acuminata_otuPS36"
#> [23] "Ceiba_crispiflora_otuPS41"
#> [24] "Pachira_aquatica_otuPS35"
#> [25] "Pachira_aquatica_otu376439"
#> [26] "Septotheca_tessmannii_otuPS11"
#> [27] "Triplochiton_zambesiacus_otuPS43"
#> [28] "Heritiera_elata_otu376445"
#> [29] "Heritiera_littoralis_otu376454"
#> [30] "Heritiera_aurea_otu376446"
#> [31] "Heritiera_simplicifolia_otu376427"
#> [32] "Heritiera_aurea_otu376435"
#> [33] "Brachychiton_acerifolius_otu376453"
#> [34] "Brachychiton_acerifolius_otu376441"
#> [35] "Acropogon_bullatus_otu376442"
#> [36] "Acropogon_dzumacensis_otu376429"
#> [37] "Franciscodendron_laurifolium_otu376443"
#> [38] "Argyrodendron_peralatum_otu376431"
#> [39] "Sterculia_balanghas_otu376450"
#> [40] "Sterculia_tragacantha_otu376428"
#> [41] "Sterculia_tragacantha_otuPS28"
#> [42] "Sterculia_stipulata_otu376440"
#> [43] "Sterculia_coccinea_otu376436"
#> [44] "Sterculia_parviflora_otu376444"
#> [45] "Hildegardia_barteri_otu376432"
#> [46] "Firmiana_malayana_otu376449"
#> [47] "Firmiana_platanifolia_otu376425"
#> [48] "Hildegardia_populifolia_otu376448"
#> [49] "Scaphium_linearicarpum_otu376434"
#> [50] "Scaphium_macropodum_otu376426"
#> [51] "Pterocymbium_tinctorium_otu376433"
#> [52] "Scaphium_macropodum_otu376451"
#> [53] "Octolobus_spectabilis_otu376447"
#> [54] "Cola_acuminata_otu376437"
#> [55] "Pterygota_alata_otu376438"
```

Trick for the cophylo titles and margins from https://cran.r-project.org/web/packages/phangorn/vignettes/IntertwiningTreesAndNetworks.html
