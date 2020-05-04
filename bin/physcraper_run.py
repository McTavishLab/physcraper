#!/usr/bin/env python
import argparse
import sys
import physcraper
from physcraper.opentree_helpers import scraper_from_opentree

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id")
parser.add_argument("-t","--tree_id", help="OpenTree tree id")
parser.add_argument("-tl", "--tree_link", help="Link to tree to update on OpenTree")
parser.add_argument("-a","--alignment", help="path to alignment")
parser.add_argument("-as","--aln_schema", help="alignment schema (nexus or fasta)")
parser.add_argument("-db", "--blast_db", help="local download of blast database")
parser.add_argument("-o","--output", help="path to output directory")
parser.add_argument("-tx","--taxonomy", help="path to taxonomy")
parser.add_argument("-c","--config_file", help="path to config file")
parser.add_argument("-tb","--treebase", action="store true", help="download alignment from treebase")
parser.add_argument("-no_est","--estimate_tree", action='store true', help="run blast search and estimate tree")


#Not yet implemented
parser.add_argument("-bl","--blast_sequence", action='store true', help="run blast search, and align but do not estimate tree")
parser.add_argument("-d","--download_data", action='store true', help="write out tree and alignment, without blasting")
parser.add_argument("-bs","--bootstrap", help="number of bootstrap reps")
parser.add_argument("-gt","--get_tree", action='store true', help="get tree from opentree")
parser.add_argument("-ga","--get_aln", action='store true', help="get alignment from opentree")
parser.add_argument("-tf", "--tree_file", help="path to your tree")
parser.add_argument("-l","linker_file", help="path to .csv linking tip labels to taxon names")


args = parser.parse_args()

assert(args.output), "Output directory (-o) is required."

if args.config_file:
    conf = physcraper.ConfigObj(args.configfile)
else:
    conf = physcraper.ConfigObj()

if args.taxonomy:
    conf.taxonomy_dir = args.taxonomy


if args.blast_db:
    conf.blast_loc = "local"
    conf.blastdb = args.blast_db
    conf.set_local()

study_id =  None
if args.tree_link:
    linkl = args.tree_link.split("/")
    assert(linkl[4]=="view")
    study_id == linkl[5]
    tree_id = linkl[-1].split("=")[1]


if args.study_id:
    study_id = args.study_id

if args.alignment:
    alnfile = args.alignment
    assert(args.aln_schema)
    aln_schema = args.aln_schema


if args.treebase:
    aln_schema = "nexus"
    alnfile = "{}/{}{}.aln".format(workdir, study_id, tree_id)
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    tre, cite = get_tree_from_study(study_id, tree_id)
    sys.stdout.write("downloading best match alignment from treebase")
    tre.write(path="{}/{}{}.tre".format(workdir, study_id, tree_id), schema="nexus")
    if not os.path.exists(alnfile):
        sys.stdout.write("downloading best match alignment from treebase")
        dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id)
        aln = get_max_match_aln(tre, dataset)
        aln.write(path=alnfile, schema = aln_schema)
    else:
        sys.stdout.write("Using alignment file found at {}.".format(alnfile))

if study_id:
    scraper = scraper_from_opentree(study_id =study_id, 
                                    tree_id = tree_id, 
                                    alnfile = alnfile, 
                                    aln_schema = aln_schema,
                                    workdir = output,
                                    configfile = conf)
    sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))

if not args.no_est:
#scraper.read_blast_wrapper()
    scraper.est_full_tree()
    scraper.data.write_labelled(label='^ot:ottTaxonName')