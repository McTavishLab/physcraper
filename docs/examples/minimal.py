""" This is a minimal example taking precooked data from the tests
"""
import os
import sys
import json
import physcraper
from physcraper import OtuJsonDict, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from physcraper import opentree_helpers
from physcraper.opentree_helpers import scraper_from_opentree

configfi = "tests/data/test.config"
workdir ="physcraper_example_minimal"
aln_fi = "tests/data/tiny_test_example/test.fas"
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
# mattype="fasta"
tre_fi= "tests/data/tiny_test_example/test.tre"
# schema_trf = "newick"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
id_to_spn = "tests/data/tiny_test_example/test_nicespl.csv"

if not os.path.exists("{}".format(workdir)):
    os.makedirs("{}".format(workdir))

conf_base = ConfigObj(configfi)

ids = IdDicts(configfi)

otu_json = OtuJsonDict(id_to_spn, ids)
with open(otu_jsonfi,"w") as outfile:
    json.dump(otu_json, outfile)

ottids = [otu_json[ite]['^ot:ottId'] for ite in otu_json]

mrca = opentree_helpers.get_mrca_ott(ottids)

# Create a 'scraper' object to get data from NCBI

data_obj = generate_ATT_from_files(alnfile=aln_fi,
                             aln_schema="fasta", #mattype
                             workdir=workdir,
                             configfile=configfi,
                             treefile=tre_fi,
                             tree_schema = "newick", #schema_trf
                             otu_json=otu_jsonfi,
                             ingroup_mrca=mrca)

data_obj.tag = "minEx"

scraper = PhyscraperScrape(data_obj, ids)

sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))

scraper._blasted = 1 # this tricks PhyScraper into not blasting anything new
sys.stdout.write("Running read_blast_wrapper()...\n")
scraper.read_blast_wrapper(blast_dir=blast_dir)
sys.stdout.write("Running write_aln()...\n")
aln_path1 = scraper.data.write_aln()
#aln_path_alt = scraper.data.write_aln(filename="already_aligned_seqs.fas")
#unaln_path = scraper.write_new_seqs(filename='unaligned.fas')

sys.stdout.write("Running align_query_seqs()...\n")
scraper.align_new_seqs()
scraper.est_full_tree()
scraper.data.write_labelled(label="^ot:ottTaxonName", filename="updated_taxon_name", norepeats=True)
scraper.data.write_labelled(label="^ncbi:taxon", filename="updated_ncbi_id", norepeats=False)


# sys.stdout.write("estimating tree...")
# scraper.est_full_tree()
# scraper.data.write_labelled(label='^ot:ottTaxonName')
