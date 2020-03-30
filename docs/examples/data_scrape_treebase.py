import json

import physcraper
from physcraper.opentree_helpers import get_tree_from_study, get_dataset_from_treebase

configfi = "docs/examples/example.config"
study_id = "ot_350"
tree_id = "Tr53297"
workdir ="scrape_ot_350_treebase"


# Read in the configuration information
conf = physcraper.ConfigObj(configfi)


#Get an existing tree from the Open Tree of life, and convert it to newick format
nexson = physcraper.opentree_helpers.get_nexson(study_id)
tre, cite = get_tree_from_study(study_id, tree_id)
tre.write(path="{}{}.tre".format(study_id, tree_id), schema="nexus")


leaves = [leaf.taxon.label for leaf in tre.leaf_node_iter()]

tax_thresh = int(0.75 * len(leaves))


#Pull down an alignment from treebase.
dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id)


aln = None
##order of data matrices is arbitratry, so we choose one that matches the tree length
for mat in dataset.char_matrices:
    if len(mat) == len(tre.taxon_namespace):
        aln = mat
        #Write it out to file, os we have the 'before' alignment
        data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                            workdir=workdir,
                                                            config_obj=conf,
                                                            study_id=study_id,
                                                            tree_id=tree_id)

# If we didn't find an alignement that is an exact match, try the 1st one

if not aln:
  for mat in dataset.char_matrices:
    # To preserve taxon labels and relationships, 
    #we will combine the alignement, tree and taxon information into a single data object
    # By using the OpenTree Phylesystem API we can get the orgininal taxon names as well as the taxon mappings
    data_obj = physcraper.generate_ATT_from_phylesystem(aln=mat,
                                                        workdir=workdir,
                                                        config_obj=conf,
                                                        study_id=study_id,
                                                        tree_id=tree_id)
    if len(data_obj.aln) >= tax_thresh:
        aln = mat
        break 



assert(len(data_obj.aln) >= tax_thresh)

data_obj.write_files(treepath="before.tre", treeschema="newick", alnpath="before.fas", alnschema="fasta")
json.dump(data_obj.otu_dict, open('{}/otu_dict.json'.format(workdir), 'wb'))




# We need to create a physcraper ids object to translate between ncbi and OpenTree identifiers.
ids = physcraper.IdDicts(conf, workdir=workdir)


# Create an 'scraper' object to get data from NCBI, align it an
scraper = physcraper.PhyscraperScrape(data_obj, ids)


sys.stdout.write("{} taxa in alignement and tree\n".format(len(scraper.data.aln)))


#scraper.read_blast_wrapper()
#scraper.est_full_tree()
#scraper.data.write_labelled(label='^ot:ottTaxonName')