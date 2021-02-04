## Analysing the Physcraper results

### Rerooting the trees

*In construction*

A correctly rooted phylogeny is needed to compare relationships between two or more phylogenies.
Rooting phylogenies can be tricky. While Physcraper places a suggested root based on the taxonomic relationships in OpenTree, this root can be unreliable, especially if taxonomy is a poor fit to true evolutionary relationships.

So whenever possible, the root must be specified by the user.



### Tree comparison with Robinson-Foulds

*In construction*

  tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o pg_55_comparison


### Relabeling the trees

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('example/docs/pg_55')
    pg55.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')
