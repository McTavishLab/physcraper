#!/usr/bin/env python
import argparse
import sys
import physcraper
from physcraper.opentree_helpers import scraper_from_opentree

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id")
parser.add_argument("-t","--tree_id", help="tree id")
parser.add_argument("-a","--alignment", help="path to alignment")
parser.add_argument("-as","--aln_schema", help="alignment schema (nexus or fasta)")
parser.add_argument("-o","--output", help="path to output directory")


args = parser.parse_args()



# Create an 'scraper' object to get data from NCBI, align it an
scraper = scraper_from_opentree(study_id =args.study_id, 
                                tree_id = args.tree_id, 
                                alnfile = args.alignment, 
                                aln_schema = args.aln_schema,
                                workdir = args.output)

sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))


#scraper.read_blast_wrapper()
scraper.est_full_tree()
scraper.data.write_labelled(label='^ot:ottTaxonName')