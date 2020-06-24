## Tutorial

### Find a study with your taxon of interest

Search on OpenTree of life using 


### Run the auto update


### Compare your new tree to existing relationships


### Reroot or relabel tree

from physcraper import treetaxon
podarc = treetaxon.generate_TreeTax_from_run('test_podarcis')
podarc.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')

