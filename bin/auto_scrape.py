#!/usr/bin/env python
import argparse
import sys
import os
import physcraper
from physcraper.opentree_helpers import get_tree_from_study, scraper_from_opentree, get_max_match_aln, count_match_tree_to_aln

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id")
parser.add_argument("-t","--tree_id", help="tree id")
parser.add_argument("-o","--output", help="path to output directory")

args = parser.parse_args()



study_id =args.study_id
tree_id = args.tree_id 
workdir = args.output

#study_id = "ot_350"
#tree_id = "Tr53296"
#tree_id2 = "Tr53297"
#workdir = "tmp_treebase"



aln_schema = "nexus"
alnfile = "{}/{}{}.aln".format(workdir, study_id, tree_id)

if not os.path.exists(workdir):
        os.makedirs(workdir)

#Get an existing tree from the Open Tree of life, and convert it to newick format
tre, cite = get_tree_from_study(study_id, tree_id)
tre.write(path="{}/{}{}.tre".format(workdir, study_id, tree_id), schema="nexus")

if not os.path.exists(alnfile):
    sys.stdout.write("downloading best match alignment from treebase")
    dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id)
    aln = get_max_match_aln(tre, dataset)
    aln.write(path=alnfile, schema = aln_schema)



# Create an 'scraper' object to get data from NCBI, align it an
scraper = scraper_from_opentree(study_id =study_id, 
                                tree_id = tree_id, 
                                alnfile = alnfile, 
                                aln_schema = aln_schema,
                                workdir = workdir)

sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))


#scraper.read_blast_wrapper()
scraper.est_full_tree()
scraper.data.write_labelled(label='^ot:ottTaxonName')