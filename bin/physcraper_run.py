#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  Physcraper Library.
##
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Luna L. Sanchez Reyes, Martha Kandziora, Emily Jane McTavish. (2020).
##     Physcraper: A Python package for continual update of evolutionary
##     estimates using the Open Tree of Life. BioRxiv 2020.09.15.299156.
##     doi: doi.org/10.1101/2020.09.15.299156
##
##############################################################################

"""
This module handles the code to generate the Physcraper command line code interface
to update a tree from an existing alignment, by:
getting a tree from an OpenTree study ID or from a local file, getting the
corresponding alignment form treeBASE or from a local file, identifying the MRCA,
running BLAST, downloading new matching sequences, aligning them and inferring a tree.
"""

import argparse
import os
import sys
import physcraper
from physcraper.opentree_helpers import (get_tree_from_study,
                                         scraper_from_opentree,
                                         get_max_match_aln,
                                         generate_ATT_from_phylesystem,
                                         bulk_tnrs_load)
from physcraper.aligntreetax import generate_ATT_from_run, generate_ATT_from_files

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id.")
parser.add_argument("-t","--tree_id", help="OpenTree tree id.")
parser.add_argument("-a","--alignment", help="Path to alignment.")
parser.add_argument("-as","--aln_schema", help="Alignment schema (nexus or fasta).")
parser.add_argument("-o","--output", help="Path to output directory.")
parser.add_argument("-c","--configfile", help="Path to config file.")
parser.add_argument("-tb","--treebase", action="store_true", help="Download alignment from treebase.")
parser.add_argument("-re","--reload_files",  help="Reload files and configuration from dir.")

parser.add_argument("-tf","--tree_file", help="A path to a tree file.")
parser.add_argument("-tfs","--tree_schema", help="Tree file format schema.")
parser.add_argument("-ti","--taxon_info", help="Taxon info file from OpenTree TNRS.")



parser.add_argument("-tag","--tag", help="Gene name or other name specifier.")
msg1 = "Taxonomic id to constrain a BLAST search. "
msg2 = "Format is ott:123 or ncbi:123. Default will use ingroup of tree on OpenTree, or MRCA of input tips."
parser.add_argument("-st","--search_taxon", help=msg1 + msg2)
msg3 = "MAXimum "
msg4 = "MINimum "
parser.add_argument("-spn","--species_number", help=msg3 + "number of sequences to include per species.")
parser.add_argument("-bl","--block_list", nargs='+', help="NCBI accession numbers to exclude.")

parser.add_argument("-tp","--trim_perc", help= msg4 + "percentage of sequences end of alignments.")

msg5 = "relative length of added sequences, compared to input sequence length."
parser.add_argument("-rlmax","--relative_length_max", help=msg3 + msg5)
parser.add_argument("-rlmin","--relative_length_min", help=msg4 + msg5)

parser.add_argument("-r","--repeat",  action='store_true', help="Repeat BLAST search until no new sequences are found.")
parser.add_argument("-db", "--blast_db", help="Local download of BLAST database.")
parser.add_argument("-u", "--blast_url", help="A URL for your own mirrored BLAST database.")
parser.add_argument("-e","--email", help="email address for NCBI BLAST searches.")
parser.add_argument("-ak","--api_key", help="NCBI API key.")

parser.add_argument("-ev","--eval", help="BLAST e-value cutoff.")
parser.add_argument("-hl","--hitlist_len", help="Number of blast searches to save per taxon.")

parser.add_argument("-nt","--num_threads", help="Number of threads to use in processing.")
parser.add_argument("-de","--delay", help="how long to wait before blasting the same sequence again")
msg6 = "Just checks the input files. Does not estimate tree, nor runs BLAST search."
parser.add_argument("-no_est","--no_estimate", action='store_true', help=msg6)

parser.add_argument("-bs","--bootstrap_reps", help="Number of bootstrap repetitions.")

parser.add_argument("-tx","--taxonomy", help="Path to taxonomy.")

parser.add_argument("-v","--verbose", action="store_true", help="Verbose mode.")



# Not yet implemented
# parser.add_argument("-tl", "--tree_link", help="Link to tree to update on OpenTree")
# msg7 = "run blast search, and align but do not estimate tree"
# parser.add_argument("-bl","--blast_sequence", action='store_true', help=msg7)
# msg8 = "write out tree and alignment, without blasting"
# parser.add_argument("-d","--download_data", action='store_true', help=msg8)
# parser.add_argument("-bs","--bootstrap", help="number of bootstrap reps")
# parser.add_argument("-gt","--get_tree", action='store_true', help="get tree from opentree")
# parser.add_argument("-ga","--get_aln", action='store_true', help="get alignment from opentree")
# parser.add_argument("-tf", "--tree_file", help="path to your tree")
# parser.add_argument("-l","--linker_file", help="path to .csv linking tip labels to taxon names")


args = parser.parse_args()


if args.verbose:
    physcraper.scrape.set_verbose()


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

try:
    assert args.output
    workdir = args.output
    if not os.path.exists(workdir):
        os.makedirs(workdir)
except AssertionError:
    sys.stderr.write("ERROR: Output directory (-o) is required.\n")
    sys.exit(-1)


#set configuration
if args.configfile:
    conf = physcraper.ConfigObj(args.configfile)
elif args.reload_files:
    files = [f for f in os.listdir(args.reload_files)]
    for file in files:
        if file.startswith('run_'):
            tag = file.split('.')[0].replace('run_', '')
    configfile = "{}/run_{}/run.config".format(args.reload_files, tag)
    conf = physcraper.ConfigObj(configfile)
    conf.workdir = args.output
else:
    conf = physcraper.ConfigObj()

if args.taxonomy:
    conf.taxonomy_dir = args.taxonomy

if args.blast_db:
    conf.blast_loc = "local"
    conf.blastdb = args.blast_db
    conf.set_local()

if args.eval:
    conf.e_value_thresh = float(args.eval)

if args.hitlist_len:
    conf.hitlist_size = int(args.hitlist_len)

if args.trim_perc:
    conf.trim_perc = float(args.trim_perc)

if args.relative_length_max:
    conf.maxlen = float(args.relative_length_max)

if args.relative_length_min:
    conf.minlen = float(args.relative_length_min)

if args.species_number:
    conf.spp_threshold = int(args.species_number)

if args.num_threads:
    conf.num_threads = int(args.num_threads)

if args.delay:
    conf.delay = int(args.delay)

if args.email:
    conf.email = args.email

if args.api_key:
    conf.api_key = args.api_key

sys.stdout.write("Configuration Settings\n")
sys.stdout.write(conf.config_str()+"\n")

study_id =  None
#if args.tree_link:
#    linkl = args.tree_link.split("/")
#    assert(linkl[5]=="view")
#    study_id = linkl[6]
#    tree_id = linkl[-1].split("=")[1]


if args.study_id or args.tree_id:
    try:
        study_id = args.study_id
        tree_id = args.tree_id
        assert study_id
        assert tree_id
    except AssertionError:
        sys.stderr.write("ERROR: To get tree from OpenTree, specify both -s [study_id] and -t [tree_id]\n")
        sys.exit(-1)


if args.alignment:
    alnfile = args.alignment
    try:
        assert args.aln_schema
        aln_schema = args.aln_schema
    except AssertionError:
        sys.stderr.write("ERROR: Specify alignment format using -as [fasta or nexus]\n")
        sys.exit(-1)


if args.treebase:
    aln_schema = "nexus"
    alnfile = "{}/{}{}.aln".format(workdir, study_id, tree_id)
    tre, cite = get_tree_from_study(study_id, tree_id)
    tre.write(path="{}/OpenTree_{}{}.tre".format(workdir, study_id, tree_id), schema="nexus")
    if not os.path.exists(alnfile):
        sys.stdout.write("downloading best match alignment from treebase to {}\n".format(alnfile))
        dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id)
        aln = get_max_match_aln(tre, dataset)
        aln.write(path=alnfile, schema = aln_schema)
    else:
        sys.stdout.write("Using alignment file found at {}.\n".format(alnfile))

# TODO Luna: add a test for this chunk of code:
search_ott_id = None
if args.search_taxon:
    ids = physcraper.IdDicts(conf)
    if args.search_taxon.lower().startswith('ott'):
        search_ott_id = args.search_taxon.split(':')[1]
    elif args.search_taxon.lower().startswith('ncbi'):
        ncbi_id = inst(args.search_taxon.split(':')[1])
        search_ott_id = ids.ncbi_ott[ncbi_id]
    else:
        sys.stderr.write("search taxon id must be in format ott:123 or ncbi:123\n")

if study_id:
    tre, cite = get_tree_from_study(study_id, tree_id)
    tre.write(path="{}/OpenTree_{}{}.tre".format(workdir, study_id, tree_id), schema="nexus")
    if search_ott_id:
        data_obj = generate_ATT_from_phylesystem(study_id =study_id,
                                                 tree_id = tree_id,
                                                 alnfile = alnfile,
                                                 aln_schema = aln_schema,
                                                 workdir = workdir,
                                                 configfile = conf,
                                                 search_taxon = search_ott_id)
        scraper = physcraper.PhyscraperScrape(data_obj, ids)
        scraper.data.write_otus(schema='table', direc=scraper.inputsdir)
        scraper.data.write_otus(schema='json', direc=scraper.inputsdir)
    else:
        scraper = scraper_from_opentree(study_id =study_id,
                                        tree_id = tree_id,
                                        alnfile = alnfile,
                                        aln_schema = aln_schema,
                                        workdir = workdir,
                                        configfile = conf)
        scraper.data.write_otus(schema='table', direc=scraper.inputsdir)
        scraper.data.write_otus(schema='json', direc=scraper.inputsdir)
    sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))


if args.tree_file:
    treefile = args.tree_file
    ast1 = "When passing in a treefile, a tree schema is required.\n"
    ast2 = "When passing in a treefile, a taxon mapping from \
            https://tree.opentreeoflife.org/curator/tnrs/ is required.\n"
    ast3 = "JSON format file required for taxon info from \
            https://tree.opentreeoflife.org/curator/tnrs/.\n"
    assert(args.tree_schema), ast1
    assert(args.taxon_info), ast2
    assert(args.taxon_info.split('.')[-1]=='json'), ast3
    if args.tag:
        tag = args.tag
    elif args.alignment:
        tag = args.alignment.split('/')[-1].split('.')[0]
    otu_dict = bulk_tnrs_load(args.taxon_info)
    if args.search_taxon:
        search_taxon = args.search_taxon
    else:
        search_taxon = None
    data_obj = generate_ATT_from_files(workdir= workdir,
                                       configfile=conf,
                                       alnfile = alnfile,
                                       aln_schema = aln_schema,
                                       treefile = treefile,
                                       otu_json = otu_dict,
                                       tree_schema = args.tree_schema,
                                       search_taxon=search_ott_id)
    ids = physcraper.IdDicts(conf)
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    scraper.data.write_otus(schema='table', direc=scraper.inputsdir)
    scraper.data.write_otus(schema='json', direc=scraper.inputsdir)

#    sys.stdout.write("Read in tree {} taxa in alignment and tree\n".format(len(scraper.data.aln)))


if args.reload_files:
    if args.tag:
        tag = args.tag
    elif args.alignment:
        tag = args.alignment.split('/')[-1].split('.')[0]
    data_obj = generate_ATT_from_run(args.reload_files, configfile=conf)
    data_obj.workdir = workdir
    ids = physcraper.IdDicts(conf)
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    sys.stdout.write("Reloaded {} taxa in alignment and tree\n".format(len(scraper.data.aln)))


if args.bootstrap_reps:
    boot_reps = args.bootstrap_reps
else:
    boot_reps = 100

if args.block_list:
    for acc in args.block_list:
        scraper.blocklist.append(acc)
        sys.stdout.write("Excluding accession numbers {}\n".format(','.join(scraper.blocklist)))


run = 0
if args.repeat:
    besttreepath = None
    to_be_blasted = ['first_pass']
    rundir_base = scraper.rundir
    while len(to_be_blasted) >= 1:
        run += 1
        scraper.run_blast_wrapper()
        besttreepath = scraper.est_full_tree()
        if besttreepath:
            prev_besttreepath = besttreepath
            scraper.replace_tre(besttreepath)
            scraper.data.write_labelled(filename="run_{}".format(run), 
                                        label='^ot:ottTaxonName', 
                                        direc=scraper.outputsdir)
            scraper.data.write_otus(schema='json', 
                                    direc=scraper.rundir)
            scraper.data.write_otus(schema='json', 
                                    direc=scraper.outputsdir)   
            scraper.data.write_files(direc=scraper.outputsdir)         
            new_rundir = "{}_run{}".format(rundir_base, 
                                           run)
            prev_rundir = scraper.rundir
            scraper.rundir = new_rundir
            os.mkdir(scraper.rundir)
            new_to_be_blasted = []
            for otu in scraper.data.aln:
                cond1 = scraper.data.otu_dict[otu.label]['^physcraper:ingroup']
                cond2 = scraper.data.otu_dict[otu.label]['^physcraper:last_blasted'] is None
                if (cond1 and cond2):
                    new_to_be_blasted.append(otu.label)
            to_be_blasted = new_to_be_blasted
        elif besttreepath is None:
          #`  os.rmdir(scraper.rundir)
            scraper.rundir = prev_rundir
            updated_alnfi = "{}/physcraper_{}.fas".format(prev_rundir, scraper.data.tag)
            bootpath = scraper.calculate_bootstrap(alignment = updated_alnfi,
                                                   num_reps = boot_reps)
            sumtreepath = scraper.summarize_boot(prev_besttreepath,
                                                 bootpath)
            scraper.replace_tre(sumtreepath,
                                schema="nexus")
            scraper.data.write_files(direc=scraper.outputsdir)
            scraper.data.write_otus(schema='table',
                                    direc=scraper.inputsdir)
            scraper.data.write_labelled(filename='updated_taxonname',
                                        label='^ot:ottTaxonName',
                                        direc = scraper.outputsdir)
            to_be_blasted =  []
        else:
            sys.stderr.write("unexpected error")
elif not args.no_estimate:
#scraper.read_blast_wrapper()
    scraper.calculate_final_tree(boot_reps = boot_reps)
    scraper.data.write_labelled(label='^ot:ottTaxonName',
                                direc=scraper.outputsdir)
    scraper.data.write_otus(schema='table', 
                            direc=scraper.outputsdir)
    scraper.data.write_otus(schema='json', 
                            direc=scraper.rundir)
    scraper.data.write_otus(schema='json', 
                            direc=scraper.outputsdir)   
    scraper.data.write_files(direc=scraper.outputsdir)   
