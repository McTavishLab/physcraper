#!/usr/bin/env python
import json
import sys
from dendropy import Tree

otu_json = "otu_dict.json"
newick = "physcraper.tre"

otu_dict = json.load(open(otu_json,"r"))
i = 0
#'query' means it is a new otu added by physcraper
for tax in otu_dict:
    taxname = otu_dict[tax].get(u'^ot:ottTaxonName')
    status = otu_dict[tax].get(u'^physcraper:status')
    if status in ('query', 'original'):
#        sys.stdout.write("{}, {}\n".format(taxname, status))
        i += 1
print i

#read tree and rename tips to if they are queries or not

tre = Tree.get(path=newick,
               schema="newick",
               preserve_underscores=True)

    
label = u'^physcraper:status'

traits = open("output.csv","w")
for taxon in tre.taxon_namespace:
     if otu_dict[taxon.label].get('^ncbi:gi', False) or otu_dict[taxon.label].get('^ncbi:accession', False):
       status = "NEW"
    else:
       status = "original"
    traits.write("{}, {}\n".format(taxon, status))

new_names = set()
for taxon in tre.taxon_namespace:
    ### here the double names of the labelled tre files are generated.
    taxname = otu_dict[taxon.label].get(u'^ot:ottTaxonName')
#    status = otu_dict[taxon.label].get(u'^physcraper:status')
    if otu_dict[taxon.label].get('^ncbi:gi', False) or otu_dict[taxon.label].get('^ncbi:accession', False):
       status = "NEW"
    else:
       status = "original"
    new_label = "{}_{}".format(taxname,status)
    i = 0
    while new_label in new_names:
        i += 1
        new_label = "{}{}".format(i, new_label)
    new_names.add(new_label)
    taxon.label = new_label




tre.write(path="out.tre",
          schema="newick",
          unquoted_underscores=True,
          suppress_edge_lengths=False)

