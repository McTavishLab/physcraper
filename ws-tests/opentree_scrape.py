import physcraper
import sys
from dendropy import DnaCharacterMatrix


# Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "tests/output/opentree"
configfi = "tests/data/remotencbi.config"

sys.stdout.write("\nTesting 'opentree scrape (1 round)'\n")
conf = physcraper.ConfigObj(configfi, interactive=False)
print "1. {}".format(conf.email)
      
    
aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                    workdir=workdir,
                                                    config_obj=conf,
                                                    study_id=study_id,
                                                    tree_id=tree_id)

ids = physcraper.IdDicts(conf, workdir=workdir)

print "3. {}".format(ids.config.email)

data_obj.prune_short()
assert len(data_obj.aln) == 9
data_obj.write_files()

scraper = physcraper.PhyscraperScrape(data_obj, ids)
scraper.est_full_tree()
