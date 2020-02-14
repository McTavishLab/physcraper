import physcraper
import dendropy
import pickle
import sys
import os
import json
from peyotl.nexson_syntax import (
    extract_tree,
    get_subtree_otus,
    extract_otu_nexson,
    PhyloSchema,
)
import dendropy

configfi = "docs/examples/example.config"
study_id = "ot_350"
tree_id = "Tr53297"
workdir ="scrape_ot_350"


# Read in the configuration information
conf = physcraper.ConfigObj(configfi)


#Get an existing tree from the Open Tree of life, and convert it to newick format
nexson = physcraper.opentree_helpers.get_nexson(study_id, 'api')
newick = extract_tree(nexson,
                      tree_id,
                      PhyloSchema('newick',
                                   output_nexml2json='1.2.1',
                                   content="tree",
                                   tip_label="ot:originalLabel"))

tre = dendropy.Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True)


#Pull down an alignment from treebase.
dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id,
                                                                phylesystem_loc='api')

aln = None
##order of data matrices is arbitratry, so we choose one that matches the tree length
for mat in dataset.char_matrices:
    if len(mat) == len(tre.taxon_namespace):
        aln = mat

# If we didn't find an alignement that is an exact match, try the 1st one
if not aln:
  aln = dataset.char_matrices[0]

#Write it out to file, os we have the 'before' alignment
aln.write(path="{}{}.aln".format(study_id, tree_id), schema="nexus")

# If we are using an alinment we already wrote  to file earlier we can use this
#aln = dendropy.DnaCharacterMatrix.get(file=open("{}{}.aln".format(study_id, tree_id)), schema="nexus", taxon_namespace=tre.taxon_namespace)

tre.write(path="{}{}.tre".format(study_id, tree_id), schema="nexus")


# To preserve taxon labels and relationships, 
#we will combine the alignement, tree and taxon information into a single data object
# By using the OpenTree Phylesystem API we can get the orgininal taxon names as well as the taxon mappings
data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                         workdir=workdir,
                                         config_obj=conf,
                                         study_id=study_id,
                                         tree_id=tree_id)

#data_obj.write_files()
#json.dump(data_obj.otu_dict, open('{}/otu_dict.json'.format(workdir), 'wb'))

sys.stdout.write("{} taxa in alignement and tree\n".format(len(data_obj.aln)))

# We need to create a physcraper ids object to translate between ncbi and OpenTree identifiers.
ids = physcraper.IdDicts(conf, workdir=workdir)


# Create an 'scraper' object to get data from NCBI, align it an
scraper = physcraper.PhyscraperScrape(data_obj, ids)

#scraper.read_blast_wrapper()
scraper.est_full_tree()
scraper.data.write_labelled(label='^ot:ottTaxonName')