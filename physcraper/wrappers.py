#!/usr/bin/env python
import pickle
import sys
import os
import subprocess
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj,  IdDicts, PhyscraperScrape
from physcraper import FilterBlast, Settings, debug  # Concat
from dendropy import DnaCharacterMatrix
from concat import Concat


def sync_ncbi(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["rsync", "av", "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz", "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp)])
    subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


def sync_ott(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["process_ott.sh", "".format(conf.ott_ncbi)])

# generates IdDicts physcrapper class
def get_ottid(configfi, cwd):
    conf = ConfigObj(configfi)
    ids = IdDicts(conf, cwd)
    return(ids)


def standard_run(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi,
                 shared_blast_folder):
    '''looks for a json file to continue run, or builds and runs
    new analysis for as long as new seqs are found'''
    debug('Debugging mode is on')

    conf = ConfigObj(configfi)
    if os.path.isfile("{}/att_checkpoint.p".format(workdir)):
        sys.stdout.write("Reloading data object from pickle file\n")
        data_obj = pickle.load(open("{}/att_checkpoint.p".format(workdir), "rb"))
#        scraper.repeat = 1
    else:
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)
        aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_phylesystem(aln=aln,
                                                 workdir=workdir,
                                                 study_id=study_id,
                                                 tree_id=tree_id,
                                                 phylesystem_loc=conf.phylesystem_loc)
        # Mapping identifiers between OpenTree and NCBI requires and identifier dict object
        ids = IdDicts(conf, workdir="example")
        # Prune sequences below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName')
        data_obj.write_otus('otu_info', schema='table')
        data_obj.dump()
        # Mapping identifiers between OpenTree and NCBI requires and identifier dict object
    if os.path.isfile(conf.id_pickle):
        sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
#       ids = jsonpickle.decode(thawed_id)
#       scraper.repeat = 1
        ids = pickle.load(open(conf.id_pickle, "rb"))
    else:
        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        ids.dump()
    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    scraper = PhyscraperScrape(data_obj, conf)
    # run the ananlyses
    # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
    # scraper.blast_subdir = "/home/martha/physcraper/phyruns/blast_runs/"
    if shared_blast_folder:
        scraper.blast_subdir = shared_blast_folder
    else:
        shared_blast_folder = None
    scraper.run_blast()
    scraper.read_blast(blast_dir= shared_blast_folder)
    scraper.remove_identical_seqs()
    scraper.generate_streamed_alignment()
    while scraper.repeat == 1:
        scraper.data.write_labelled(label='^ot:ottTaxonName')
        scraper.data.write_otus("otu_info", schema='table')
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.run_blast()
        scraper.read_blast(blast_dir= shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    return


def own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 sp_info_jsonfi,
                 configfi,
                 ingroup_mrca = None,
		         shared_blast_folder=None):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found'''

    debug('Debugging mode is on')

    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: ATT\n")
        scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
        scraper.repeat = 1    
    else:   
        # sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)

        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                             mattype=mattype,
                                             workdir=workdir,
                                             treefile=trfn,
                                             schema_trf = schema_trf,
                                             otu_json=sp_info_jsonfi,
                                             ingroup_mrca=None)

        # Prune sequnces below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled( label='^ot:ottTaxonName')
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        scraper =  PhyscraperScrape(data_obj, ids)
        print(scraper.data.aln.taxon_namespace)
        print(scraper.data.tre.taxon_namespace)

        # run the analyses
        # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
        # scraper.blast_subdir = "/home/martha/physcraper/phyruns/blast_runs/"
        # scraper.blast_subdir = "/home/mkandziora/shared_runs/"
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.run_blast()
        scraper.read_blast(blast_dir= shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    while scraper.repeat == 1: 
        scraper.run_blast()
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.read_blast(blast_dir= shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    return 1 


def filter_OTOL(study_id,
                 tree_id,
                 seqaln,
                 workdir,
                 configfi,
                 treshold,
                 selectby,
                 downtorank,
                 blacklist,
                 add_local_seq,
                 id_to_spn_addseq_json):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.'''
    debug('Debugging mode is on')

    # if _DEBUG_MK == 1:
    #     random.seed(3269235691)
    print(workdir)
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        filteredScrape.repeat = 1   
    else:   
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        # print("config")
        debug(dir(conf))
        debug(conf.email)

        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_phylesystem(seqaln,
                                  workdir,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api')

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled( label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()

        ids = IdDicts(conf, workdir=workdir)

        #Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape =  FilterBlast(data_obj, ids)


        filteredScrape.blacklist = blacklist

# filteredScrape.write_otu_info(downtorank)
        if add_local_seq is not None:
            debug("will add local sequences now")
            filteredScrape.add_local_seq(add_local_seq, id_to_spn_addseq_json)
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
        #run the ananlyses
        sys.stdout.write("BLASTing input sequences\n")
        # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
        #filteredScrape.blast_subdir = "/home/martha/physcraper/phyruns/blast_runs/"

        filteredScrape.run_blast()
        filteredScrape.read_blast(blast_dir= shared_blast_folder)
        sys.stdout.write("remove idential sequences\n")
        filteredScrape.remove_identical_seqs()
        filteredScrape.dump()
        debug(treshold)
        sys.stdout.write("Filter the sequences\n")
        if treshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        debug("from replace to streamed aln")
        sys.stdout.write("calculate the phylogeny\n")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
    while filteredScrape.repeat == 1:
        # number_rounds += 1
        filteredScrape.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        filteredScrape.data.write_otus("otu_info", schema='table')
        sys.stdout.write("BLASTing input sequences\n")
        filteredScrape.run_blast()
        filteredScrape.read_blast(blast_dir= shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        sys.stdout.write("Filter the sequences\n")
        if treshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.data.reconcile(seq_len_perc=0.75)
        sys.stdout.write("calculate the phylogeny\n")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        filteredScrape.write_otu_info(downtorank)
        return filteredScrape



def filter_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downtorank,
                 spInfoDict,
                 blacklist,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi,
                 ingroup_mrca=None,
                 shared_blast_folder=None):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.'''
    debug('Debugging mode is on')

    # debug(shared_blast_folder)
    # debug(some)
    # # if _DEBUG_MK == 1:
    #     random.seed(3269235691)
    print(workdir)
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        filteredScrape.repeat = 1   
    else:   
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        # print("config")
        debug(dir(conf))
        debug(conf.email)

        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             treefile=trfn,
                             schema_trf=schema_trf,
                             otu_json=spInfoDict,
                             #email=conf.email,
                             ingroup_mrca=None)

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled( label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()

        ids = IdDicts(conf, workdir=workdir)

        #Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape =  FilterBlast(data_obj, ids)
        filteredScrape.blacklist = blacklist

        # filteredScrape.write_otu_info(downtorank)
        debug(add_local_seq)
        debug(add_local_seq != None)
        if add_local_seq is not None:
            filteredScrape.localblast = True
            debug(filteredScrape.localblast)
        if filteredScrape.localblast == True:

            debug("will add local sequences now")
            filteredScrape.localblast = True
            filteredScrape.write_unpl_lblastdb(add_local_seq)
            #filteredScrape.run_local_blast()
            filteredScrape.run_blast()
            print(id_to_spn_addseq_json)
            filteredScrape.local_otu_json = id_to_spn_addseq_json

            filteredScrape.read_blast()
            sys.stdout.write("remove idential sequences\n")
            filteredScrape.remove_identical_seqs()

            filteredScrape.generate_streamed_alignment()
            # print(some)
            filteredScrape.localblast = False

        else:
            #run the ananlyses
            sys.stdout.write("BLASTing input sequences\n")
            # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
            print(shared_blast_folder)
            if shared_blast_folder:
                filteredScrape.blast_subdir = shared_blast_folder
            else:
                shared_blast_folder = None
            print(shared_blast_folder)

            filteredScrape.run_blast()
            filteredScrape.read_blast(blast_dir= shared_blast_folder)
            sys.stdout.write("remove idential sequences\n")
            filteredScrape.remove_identical_seqs()
            filteredScrape.dump()
            debug(treshold)
            sys.stdout.write("Filter the sequences\n")
            if treshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            debug("from replace to streamed aln")
            sys.stdout.write("calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
        while filteredScrape.repeat == 1:
            # number_rounds += 1
            filteredScrape.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
            filteredScrape.data.write_otus("otu_info", schema='table')
            sys.stdout.write("BLASTing input sequences\n")
            if shared_blast_folder:
                filteredScrape.blast_subdir = shared_blast_folder
            else:
                shared_blast_folder = None
            filteredScrape.run_blast()
            filteredScrape.read_blast(blast_dir= shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            sys.stdout.write("Filter the sequences\n")
            if treshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            filteredScrape.data.reconcile(seq_len_perc=0.75)
            sys.stdout.write("calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
            filteredScrape.write_otu_info(downtorank)
            return filteredScrape



#######################3
def make_settings_class(seqaln, mattype, trfn, schema_trf, workdir, 
                        treshold=None, selectby=None, downtorank=None, spInfoDict=None, add_local_seq=None, 
                        id_to_spn_addseq_json=None, configfi=None, blacklist=None):
    """all the settings are set here and can then be fed to the FilterClass
    """
    settings = Settings(seqaln=seqaln, mattype=mattype, trfn=trfn, schema_trf=schema_trf, workdir=workdir, 
                        treshold=treshold, selectby=selectby, downtorank=downtorank, spInfoDict=spInfoDict, 
                        add_local_seq=add_local_seq, id_to_spn_addseq_json=id_to_spn_addseq_json, configfi=configfi, blacklist=blacklist)
    return settings
    

def run_with_settings(settings):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.'''
    debug('Debugging mode is on')

    # if _DEBUG_MK == 1:
    #     random.seed(3269235691)

    if os.path.isfile("{}/scrape_checkpoint.p".format(settings.workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(settings.workdir),'rb'))
        filteredScrape.repeat = 1   
    else: 
        conf = ConfigObj(settings.configfi)
        # print("config")
        debug(dir(conf))
        debug(conf.email)

        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=settings.seqaln, 
                             mattype=settings.mattype, 
                             workdir=settings.workdir,
                             treefile=settings.trfn,
                             schema_trf=settings.schema_trf,
                             otu_json=settings.spInfoDict,
                             #email=conf.email,
                             ingroup_mrca=None)

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()

        data_obj.write_labelled( label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()
        
        ids = IdDicts(conf, workdir=settings.workdir)

        filteredScrape =  FilterBlast(data_obj, ids, settings)
        filteredScrape.write_otu_info(settings.downtorank)
        if settings.add_local_seq is not None:
            debug("will add local sequences now")
            filteredScrape.add_local_seq(settings.add_local_seq, settings.id_to_spn_addseq_json)
            # scraper.replace_new_seq()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
        #run the ananlyses
        # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
        # filteredScrape.blast_subdir = "/home/martha/physcraper/phyruns/blast_runs/"

        filteredScrape.run_blast()
        filteredScrape.read_blast(blast_dir= shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        filteredScrape.dump()
        if settings.treshold is not None:
            filteredScrape.sp_dict(settings.downtorank)
            filteredScrape.make_sp_seq_dict(treshold=settings.treshold, selectby=settings.selectby)
            filteredScrape.how_many_sp_to_keep(treshold=settings.treshold, selectby=settings.selectby)
            filteredScrape.replace_new_seq()
        debug("from replace to streamed aln")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
    while filteredScrape.repeat == 1:
        # number_rounds += 1
        filteredScrape.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        filteredScrape.data.write_otus("otu_info", schema='table')
        filteredScrape.run_blast()
        filteredScrape.read_blast(blast_dir= shared_blast_folder)
        filteredScrape.remove_identical_seqs()

        debug("make sp_dict")    
        if settings.treshold is not None:
            filteredScrape.sp_dict(settings.downtorank)
            filteredScrape.make_sp_seq_dict(treshold=settings.treshold, selectby=settings.selectby)
            filteredScrape.how_many_sp_to_keep(treshold=settings.treshold, selectby=settings.selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        filteredScrape.write_otu_info(settings.downtorank)
        return filteredScrape


def concat(genelistdict, workdir_comb, email, percentage=0.37, user_concat_fn=None):
    """This is to concatenate different physcraper runs into a single alignment and tree.
    genelistdict is a dict with gene names as key and the corresponding workdir
    """
   
    concat = Concat(workdir_comb, email)
    concat.concatfile = user_concat_fn
    # print(genelistdict)
    for item in genelistdict.keys():
        concat.load_single_genes(genelistdict[item]["workdir"], genelistdict[item]["pickle"], item)

    concat.combine()
    concat.sp_seq_counter()
    sp_to_keep = concat.sp_to_keep()
    concat.get_largest_tre()
    concat.make_sp_gene_dict(sp_to_keep)
    concat.make_alns_dict()
    concat.concatenate_alns()
    concat.get_short_seq_from_concat(percentage)
    concat.remove_short_seq()
    concat.make_concat_table()
    concat.write_partition()
    concat.place_new_seqs()
    concat.est_full_tree()
    concat.calculate_bootstrap()
    return concat
    