## Analysing the Physcraper results

### Tree comparison with Robinson-Foulds

*In construction*

### Relabeling the trees

    from physcraper import treetaxon
    pg55 = treetaxon.generate_TreeTax_from_run('example/docs/pg_55')
    pg55.write_labelled(label='^ot:ottTaxonName', norepeats=False, path='test_podarcis/repeats.tre')

### Rerooting the trees

*In construction*
