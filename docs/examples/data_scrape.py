import sys
import physcraper
from physcraper.opentree_helpers import scraper_from_opentree

configfi = "docs/examples/example.config"
study_id = "ot_350"
tree_id = "Tr53297"
workdir ="physcraper_example_ot_350"
aln_fi = "docs/examples/{}{}.aln".format(study_id, tree_id)


# Create an 'scraper' object to get data from NCBI, align it an
scraper = scraper_from_opentree(study_id = study_id,
                                tree_id = tree_id,
                                alnfile = aln_fi,
                                aln_schema = "nexus",
                                configfile = configfi,
                                workdir = workdir)

sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))


#scraper.read_blast_wrapper()
scraper.est_full_tree()
scraper.data.write_labelled(label='^ot:ottTaxonName')
