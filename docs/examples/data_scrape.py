import sys
import physcraper
from physcraper.opentree_helpers import scraper_from_opentree

configfi = "docs/examples/example.config"
study_id = "ot_350"
tree_id = "Tr53297"
workdir ="scrape_ot_350_compact1"


aln_fi = "docs/examples/{}{}.aln".format(study_id, tree_id)

#data_obj.write_files()
#json.dump(data_obj.otu_dict, open('{}/otu_dict.json'.format(workdir), 'wb'))


# Create an 'scraper' object to get data from NCBI, align it an
scraper = scraper_from_opentree(study_id = study_id, tree_id = tree_id, alnfile = aln_fi, aln_schema = "nexus", configfile = configfi, workdir = workdir)

sys.stdout.write("{} taxa in alignement and tree\n".format(len(scraper.data.aln)))


#scraper.read_blast_wrapper()
scraper.est_full_tree()
scraper.data.write_labelled(label='^ot:ottTaxonName')