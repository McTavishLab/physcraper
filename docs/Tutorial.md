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


The script `physcraper_run.py` wraps together linking the trees an alignment, blasting, aligning sequences, and 
inferring an updated tree.  
The blast search part of updating trees take a long time! (this analysis took around 12 hours)


    physcraper_run.py -s pg_55 -t tree5864 -tb -r -o pg_55


We have put example outputs from this command in docs/examples/pg_55

### Output stucture

The analysis folder has several subdirectories.
each folder is labelled with a 'tag', which by default is the alignment name, but can be set in the physcraper_run.py arguments.

The stucture constsins of 

    -  inputs_tag
        This conatins the input tree and alignemnt, and the mapping of teh labels to taxa saved as otu_info.csv
    -  blast_run_tag
    This directory contains the blast results for each tip in the tree
    -  run_tag
    This is where intermediate processing files are stored
    - outputs_tag
    The final tree and alignment are saved here

the final inferences will be in 

output_directory/outputs_tag

The key files are:

The inferred tree, physcraper.tre and teh updated alignment, physcraper.tre.

The tips on the tree and the alignemnt are labeled with unique identifiers as labels.

The data associacted with each of these labels is described in human readable format in out_info.csv. 
The canonical otu info is also stored in json format in 


This same tree and alignment are also output with teh taxon names as labels, saves in outputs as updated_taxonname.tre and .fas.


To further explore the data associated with the tips on this tree see DataExploration.md 



### Compare your new tree to existing relationships

    tree_comparison.py -d docs/examples/pg_55/ -og otu376420 otu376439 otu376452 -o pg_55_comparison



### Reroot or relabel tree

    from physcraper import treetaxon
    podarc = treetaxon.generate_TreeTax_from_run('test_podarcis')
    podarc.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')



##Example with Data Dryad chiroptera gene trees???
