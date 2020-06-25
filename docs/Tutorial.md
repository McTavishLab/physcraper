## Tutorial

### Find a study with your taxon of interest

Search on OpenTree of life using your taxon of interest, e.g. 'Malvaceae'

    find_trees.py --taxon_name "Malvaceae"

Pick a study!

We chose Wilkie et al (2006) https://tree.opentreeoflife.org/curator/study/view/pg_55

While this study was focussed on the family "Sterculiacea", 
phylogenetic inference have suggested that this family is not monophyletic.

https://tree.opentreeoflife.org/opentree/argus/ottol@996482

In order to further assess, lets update the tree!


### Run the auto update

The blast search part of updating trees take a long time! (hours)


    physcraper_run.py -s pg_55 -t tree5864 -tb -r -o pg_55


We have put example outputs from this command in docs/examples/pg_55

### Compare your new tree to existing relationships

    tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o tmp_comp



### Reroot or relabel tree

    from physcraper import treetaxon
    podarc = treetaxon.generate_TreeTax_from_run('test_podarcis')
    podarc.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')

