#!/usr/bin/env python
import argparse
import os
import sys
import physcraper
from physcraper.opentree_helpers import get_tree_from_study, scraper_from_opentree, get_max_match_aln, count_match_tree_to_aln, generate_ATT_from_phylesystem
from physcraper.aligntreetax import generate_ATT_from_run

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id")
parser.add_argument("-t","--tree_id", help="OpenTree tree id")
parser.add_argument("-tl", "--tree_link", help="Link to tree to update on OpenTree")
parser.add_argument("-a","--alignment", help="path to alignment")
parser.add_argument("-as","--aln_schema", help="alignment schema (nexus or fasta)")
parser.add_argument("-db", "--blast_db", help="local download of blast database")
parser.add_argument("-o","--output", help="path to output directory")
parser.add_argument("-tx","--taxonomy", help="path to taxonomy")
parser.add_argument("-c","--configfile", help="path to config file")
parser.add_argument("-e","--email", help="email address for ncbi balst searches")
parser.add_argument("-re","--reload_files",  help="reload files and configureation from dir")
parser.add_argument("-tag","--tag", help="gene name or other specifier")
parser.add_argument("-tb","--treebase", action="store_true", help="download alignment from treebase")
parser.add_argument("-no_est","--no_estimate_tree", action='store_true', help="run blast search and estimate tree")
parser.add_argument("-ev","--eval", help="blast evalue cutoff")
parser.add_argument("-hl","--hitlist_len", help="number of blast searches to save per taxon")
parser.add_argument("-tp","--trim_perc", help="minimum percentage of seqs end of alignemnts")
parser.add_argument("-rl","--relative_length", help="max relative length of added seqs, compared to input aln len")
parser.add_argument("-spn","--species_number", help="max number of seqs to include per species")
parser.add_argument("-nt","--num_threads", help="number of threads to use in processing")
parser.add_argument("-de","--delay", help="how long to wait before blasting the same sequence again")
parser.add_argument("-st","--search_taxon", help="taxonomic id to constrain blast search. format ott:123 or ncbi:123. Deafult will use ingroup of tree on OpenTree, or MRCA of input tips ")





#Not yet implemented
#parser.add_argument("-bl","--blast_sequence", action='store_true', help="run blast search, and align but do not estimate tree")
#parser.add_argument("-d","--download_data", action='store_true', help="write out tree and alignment, without blasting")
#parser.add_argument("-bs","--bootstrap", help="number of bootstrap reps")
#parser.add_argument("-gt","--get_tree", action='store_true', help="get tree from opentree")
#parser.add_argument("-ga","--get_aln", action='store_true', help="get alignment from opentree")
#parser.add_argument("-tf", "--tree_file", help="path to your tree")
#parser.add_argument("-l","--linker_file", help="path to .csv linking tip labels to taxon names")


args = parser.parse_args()

assert(args.output), "Output directory (-o) is required."
workdir = args.output

#set configuration
if args.configfile:
    conf = physcraper.ConfigObj(args.configfile)
elif args.reload_files:
    configfile = "{}/run.config".format(args.reload_files)
    conf = physcraper.ConfigObj(configfile)
    sys.stdout.write("Using config file {}\n".format(configfile))
else:
    conf = physcraper.ConfigObj()

if args.taxonomy:
    conf.taxonomy_dir = args.taxonomy

if args.blast_db:
    conf.blast_loc = "local"
    conf.blastdb = args.blast_db
    conf.set_local()

if args.eval:
    conf.e_value_thresh = args.eval

if args.hitlist_len:
    conf.hitlist_size = args.hitlist_len

if args.trim_perc:
    conf.trim_perc = args.trim_perc

if args.relative_length:
    conf.maxlen = args.relative_length

if args.species_number:
    conf.spp_threshold = int(args.species_number)

if args.num_threads:
   conf.num_threads = args.num_threads

if args.delay:
    conf.delay = args.delay

if args.email:
    conf.email = args.email

sys.stdout.write("Configuration Settings\n")
sys.stdout.write(conf.config_str())

study_id =  None
if args.tree_link:
    linkl = args.tree_link.split("/")
    assert(linkl[4]=="view")
    study_id == linkl[5]
    tree_id = linkl[-1].split("=")[1]

if args.study_id:
    study_id = args.study_id
    tree_id = args.tree_id

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
    tre.write(path="{}/{}{}.tre".format(workdir, study_id, tree_id), schema="nexus")
    if not os.path.exists(alnfile):
        sys.stdout.write("downloading best match alignment from treebase\n")
        dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id)
        aln = get_max_match_aln(tre, dataset)
        aln.write(path=alnfile, schema = aln_schema)
    else:
        sys.stdout.write("Using alignment file found at {}.\n".format(alnfile))

if study_id:
    if args.search_taxon:
        ids = physcraper.IdDicts(conf)
        if args.search_taxon.startswith('ott'):
            ott_id = args.search_taxon.split(':')[1]
        elif args.search_taxon.startswith('ncbi'):
            ncbi_id = inst(args.search_taxon.split(':')[1])
            ott_id = ids.ncbi_ott[ncbi_id]
        else:
            sys.stderr.write("search taxon id must be in format ott:123 or ncbi:123\n")
        data_obj = generate_ATT_from_phylesystem(study_id =study_id, 
                                                tree_id = tree_id, 
                                                alnfile = alnfile, 
                                                aln_schema = aln_schema,
                                                workdir = workdir,
                                                configfile = conf,
                                                ingroup_mrca = ott_id)
        scraper = physcraper.PhyscraperScrape(data_obj, ids)
    else:
        scraper = scraper_from_opentree(study_id =study_id, 
                                        tree_id = tree_id, 
                                        alnfile = alnfile, 
                                        aln_schema = aln_schema,
                                        workdir = workdir,
                                        configfile = conf)
    sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))
    scraper.data.write_files()
    scraper.data.write_labelled(filename="before_labelled_{}".format(scraper.data.tag), label='^ot:ottTaxonName')
    scraper.data.write_otus(schema="json")

if args.reload_files:
    if args.tag:
        tag = args.tag
    elif args.alignment:
        tag = args.alignment.split('/')[-1].split('.')[0]
    data_obj = generate_ATT_from_run(args.reload_files, configfile=conf)
    ids = physcraper.IdDicts(conf)
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    sys.stdout.write("Reloaded {} taxa in alignment and tree\n".format(len(scraper.data.aln)))


if not args.no_estimate_tree:
#scraper.read_blast_wrapper()
    scraper.est_full_tree()
    scraper.data.write_labelled(label='^ot:ottTaxonName')