#!/usr/bin/env python
import argparse
import os
import sys
import physcraper
from physcraper.opentree_helpers import get_tree_from_study, scraper_from_opentree, get_max_match_aln, count_match_tree_to_aln, generate_ATT_from_phylesystem, bulk_tnrs_load
from physcraper.aligntreetax import generate_ATT_from_run, generate_ATT_from_files

parser = argparse.ArgumentParser()
parser.add_argument("-s","--study_id", help="OpenTree study id")
parser.add_argument("-t","--tree_id", help="OpenTree tree id")
parser.add_argument("-a","--alignment", help="path to alignment")
parser.add_argument("-as","--aln_schema", help="alignment schema (nexus or fasta)")
parser.add_argument("-o","--output", help="path to output directory")
parser.add_argument("-c","--configfile", help="path to config file")
parser.add_argument("-tb","--treebase", action="store_true", help="download alignment from treebase")
parser.add_argument("-re","--reload_files",  help="reload files and configuration from dir")

parser.add_argument("-tf","--tree_file", help="a path to a tree file")
parser.add_argument("-tfs","--tree_schema", help="tree file format schema")
parser.add_argument("-ti","--taxon_info", help="taxon info file from OpenTree TNRS")



parser.add_argument("-tag","--tag", help="gene name or other specifier")
parser.add_argument("-st","--search_taxon", help="taxonomic id to constrain blast search. format ott:123 or ncbi:123. Deafult will use ingroup of tree on OpenTree, or MRCA of input tips ")

parser.add_argument("-spn","--species_number", help="max number of seqs to include per species")
parser.add_argument("-bl","--block_list", nargs='+', help="ncbi accession numbers to exclude")

parser.add_argument("-tp","--trim_perc", help="minimum percentage of seqs end of alignemnts")
parser.add_argument("-rlmax","--relative_length_max", help="max relative length of added seqs, compared to input seq len")
parser.add_argument("-rlmin","--relative_length_min", help="min relative length of added seqs, compared to input seq len")

parser.add_argument("-r","--repeat",  action='store_true', help="repeat search until no no sequences are found")
parser.add_argument("-db", "--blast_db", help="local download of blast database")
parser.add_argument("-u", "--blast_url", help="a URL for your own mirrored blast database")
parser.add_argument("-e","--email", help="email address for ncbi blast searches")
parser.add_argument("-ak","--api_key", help="ncbi api key")

parser.add_argument("-ev","--eval", help="blast evalue cutoff")
parser.add_argument("-hl","--hitlist_len", help="number of blast searches to save per taxon")

parser.add_argument("-nt","--num_threads", help="number of threads to use in processing")
parser.add_argument("-de","--delay", help="how long to wait before blasting the same sequence again")

parser.add_argument("-no_est","--no_estimate_tree", action='store_true', help="don't estimate tree")

parser.add_argument("-bs","--bootstrap_reps", help="number of bootstrap reps")

parser.add_argument("-tx","--taxonomy", help="path to taxonomy")

parser.add_argument("-v","--verbose", action="store_true", help="OpenTree study id")



#Not yet implemented
#parser.add_argument("-tl", "--tree_link", help="Link to tree to update on OpenTree")

#parser.add_argument("-bl","--blast_sequence", action='store_true', help="run blast search, and align but do not estimate tree")
#parser.add_argument("-d","--download_data", action='store_true', help="write out tree and alignment, without blasting")
#parser.add_argument("-bs","--bootstrap", help="number of bootstrap reps")
#parser.add_argument("-gt","--get_tree", action='store_true', help="get tree from opentree")
#parser.add_argument("-ga","--get_aln", action='store_true', help="get alignment from opentree")
#parser.add_argument("-tf", "--tree_file", help="path to your tree")
#parser.add_argument("-l","--linker_file", help="path to .csv linking tip labels to taxon names")


args = parser.parse_args()


if args.verbose:
    physcraper.scrape.set_verbose()


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

try:
    assert(args.output)
    workdir = args.output
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
    conf.e_value_thresh = args.eval

if args.hitlist_len:
    conf.hitlist_size = args.hitlist_len

if args.trim_perc:
    conf.trim_perc = args.trim_perc

if args.relative_length_max:
    conf.maxlen = args.relative_length_max

if args.relative_length_min:
    conf.minlen = args.relative_length_min

if args.species_number:
    conf.spp_threshold = int(args.species_number)

if args.num_threads:
   conf.num_threads = args.num_threads

if args.delay:
    conf.delay = args.delay

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
        assert(study_id)
        assert(tree_id)
    except AssertionError:
        sys.stderr.write("ERROR: To get tree from OpenTree, specify both -s [study_id] and -t [tree_id]\n")
        sys.exit(-1)


if args.alignment:
    alnfile = args.alignment
    try:
        assert(args.aln_schema)
        aln_schema = args.aln_schema
    except AssertionError:
        sys.stderr.write("ERROR: Specify alignment format using -as [fasta or nexus]\n")
        sys.exit(-1)


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
                                                search_taxon = ott_id)
        scraper = physcraper.PhyscraperScrape(data_obj, ids)
    else:
        scraper = scraper_from_opentree(study_id =study_id, 
                                        tree_id = tree_id, 
                                        alnfile = alnfile, 
                                        aln_schema = aln_schema,
                                        workdir = workdir,
                                        configfile = conf)
    sys.stdout.write("{} taxa in alignment and tree\n".format(len(scraper.data.aln)))


if args.tree_file:
    treefile = args.tree_file
    assert(args.tree_schema), "When passing in a treefile, a tree schema is required\n"
    assert(args.taxon_info), "When passing in a treefile, a taxon mapping from from https://tree.opentreeoflife.org/curator/tnrs/ is required\n"
    assert(args.taxon_info.split('.')[-1]=='json'), "JSON format file required for taxon info from https://tree.opentreeoflife.org/curator/tnrs/\n"
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
                                        search_taxon=search_taxon)
    ids = physcraper.IdDicts(conf)
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
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


run = 1
if args.repeat:
    to_be_blasted = ['first_pass']
    while len(to_be_blasted) >= 1:
        run += 1
        print(scraper.blast_subdir)
        scraper.run_blast_wrapper()
        besttreepath = scraper.est_full_tree()
        scraper.replace_tre(besttreepath)
        scraper.data.write_labelled(filename="run_{}".format(run), label='^ot:ottTaxonName', direc=scraper.outputsdir)        
        scraper.data.write_otus(schema='table', direc=scraper.inputsdir)
        scraper.data.write_otus(schema='json', direc=scraper.rundir)  
        to_be_blasted = [otu.label for otu in scraper.data.aln if ((scraper.data.otu_dict[otu.label]['^physcraper:ingroup'] == True) and (scraper.data.otu_dict[otu.label]['^physcraper:last_blasted']==None))]
    scraper.calculate_final_tree(boot_reps = boot_reps)
elif not args.no_estimate_tree:
#scraper.read_blast_wrapper()
    scraper.calculate_final_tree(boot_reps = boot_reps)
    scraper.data.write_labelled(label='^ot:ottTaxonName')
    