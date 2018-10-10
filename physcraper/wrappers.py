#!/usr/bin/env python
import pickle
import sys
import os
import subprocess
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj,  IdDicts, PhyscraperScrape
from physcraper import FilterBlast, Settings, debug  # Concat
from dendropy import DnaCharacterMatrix
from concat import Concat

print("Current Version number: 09142018.0")

# TODO: we never do anything with the function nor the file
def sync_ncbi(configfi): 
    conf = ConfigObj(configfi)
    subprocess.call(["rsync", "av", "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz", "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp)])
    subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


# TODO: not used, process_ott.sh does not exist
def sync_ott(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["process_ott.sh", "".format(conf.ott_ncbi)])


# TODO: not used
# generates IdDicts physcrapper class
def get_ottid(configfi, cwd):
    conf = ConfigObj(configfi)
    ids = IdDicts(conf, cwd)
    return ids


def standard_run(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi,
                 ingroup_mrca=None,
                 shared_blast_folder=None):
    """looks for a json file to continue run, or builds and runs
    new analysis for as long as new seqs are found

    This is the wrapper function to start a PhyScraper run with tree and alignment ids from Open Tree of Life.
    You need:
         seqaln = ID of alignment file
         mattype = the format name of you alignment
         trfn = Id of phylogeny to update
         workdir = define where your analysis files shall be stored
         configfi = path to your config file
         ingroup_mrca = define the mrca, by supplying the Open Tree of Life identifier of the clade of interest

         shared_blast_folder = not necessary, if you want to share blast searches across runs (see documentation),
                                give the path to the folder with the shared runs.
    """
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
                                                 phylesystem_loc=conf.phylesystem_loc,
                                                 ingroup_mrca=ingroup_mrca)
        # Mapping identifiers between OpenTree and NCBI requires and identifier dict object
        # ids = IdDicts(conf, workdir="example")
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
    scraper = PhyscraperScrape(data_obj, ids)
    # run the analyses
    if shared_blast_folder:
        scraper.blast_subdir = shared_blast_folder
    else:
        shared_blast_folder = None
    scraper.run_blast()
    scraper.read_blast(blast_dir=shared_blast_folder)
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
        scraper.read_blast(blast_dir=shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    # scraper.write_otu_info()

    return scraper


def own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 sp_info_jsonfi,
                 configfi,
                 ingroup_mrca=None,
                 shared_blast_folder=None):
    """This is the wrapper function to start a PhyScraper run with your own data.
    You need:
         seqaln = path to sequence alignment file
         mattype = the format name of you alignment
         trfn = path to file with the phylogeny to update
         schema_trf = format type of your phylogeny
         workdir = define where your analysis files shall be stored
         sp_info_jsonfi = a json file which has the otu_dict stored, which is generated by the OtuJsonDict function
                            (usually, just leave it like it is in the example scripts.).
         configfi = path to your config file
         ingroup_mrca = not necessary, if you want to limit your run to a certain clade, give the OpenTree ID here,
                        can be obtained bu running: python scripts/get_ott.py ingroup_name
         shared_blast_folder = not necessary, if you want to share blast searches across runs (see documentation),
                                give the path to the folder with the shared runs.
    """

    debug('Debugging mode is on')

    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: ATT\n")
        scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
        scraper.repeat = 1    
    else:   
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)

        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                           mattype=mattype,
                                           workdir=workdir,
                                           treefile=trfn,
                                           schema_trf=schema_trf,
                                           otu_json=sp_info_jsonfi,
                                           ingroup_mrca=ingroup_mrca)

        # Prune sequences below a certain length threshold
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName')
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up ID dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        scraper = PhyscraperScrape(data_obj, ids)
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        # run the analyses
        scraper.run_blast()
        scraper.read_blast(blast_dir=shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    while scraper.repeat == 1: 
        scraper.run_blast()
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.read_blast(blast_dir=shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
    return 1 


def filter_OTOL(study_id,
                tree_id,
                seqaln,
                workdir,
                configfi,
                threshold,
                selectby="blast",
                downtorank=None,
                blacklist=None,
                add_unpubl_seq=None,  # path to local seq
                id_to_spn_addseq_json=None,
                ingroup_mrca=None,
                shared_blast_folder=None):
    """looks for pickeled file to continue run, or builds and runs
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output."""
    debug('Debugging mode is on')

    # if _DEBUG_MK == 1:
    #     random.seed(3269235691)
    print(workdir)
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
        filteredScrape.repeat = 1   
    else:   
        # sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)
        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_phylesystem(seqaln,
                                                 workdir,
                                                 study_id,
                                                 tree_id,
                                                 phylesystem_loc='api',
                                                 ingroup_mrca=ingroup_mrca)

        # Prune sequnces below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()

        ids = IdDicts(conf, workdir=workdir)

        # Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape = FilterBlast(data_obj, ids)
        filteredScrape.blacklist = blacklist
        if add_unpubl_seq is not None:
            filteredScrape.unpublished = True
            # debug(filteredScrape.unpublished)
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.unpublished = True
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            # filteredScrape.make_otu_dict_entry_unpubl()
            filteredScrape.run_blast()
            debug(id_to_spn_addseq_json)
            filteredScrape.data.local_otu_json = id_to_spn_addseq_json
            filteredScrape.read_blast()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False
        else:
            sys.stdout.write("BLASTing input sequences\n")
            filteredScrape.run_blast()
            filteredScrape.read_blast(blast_dir=shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            filteredScrape.dump()
            debug(threshold)
            if threshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            debug("from replace to streamed aln")
            sys.stdout.write("calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
    while filteredScrape.repeat == 1:
        filteredScrape.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        filteredScrape.data.write_otus("otu_info", schema='table')
        sys.stdout.write("BLASTing input sequences\n")
        filteredScrape.run_blast()
        filteredScrape.read_blast(blast_dir=shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        sys.stdout.write("Filter the sequences\n")
        if threshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.data.prune_short(0.75)
        sys.stdout.write("calculate the phylogeny\n")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        filteredScrape.write_otu_info(downtorank)
        return filteredScrape


def add_unpubl_to_backbone(seqaln,
                        mattype,
                        trfn,
                        schema_trf,
                        workdir,
                        threshold,
                        spInfoDict,
                        configfi,
                        selectby="blast",
                        downtorank="species",
                        blacklist=None,
                        add_unpubl_seq=None,
                        id_to_spn_addseq_json=None,
                        ingroup_mrca=None,
                        shared_blast_folder=None):
    """looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.
    It adds unpublished data to an input tree (evalue should be higher than usual).
    Backbone will not be updated
    """
    debug('Debugging mode is on')
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
        filteredScrape.repeat = 1   
    else:   
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)
       
        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                           mattype=mattype,
                                           workdir=workdir,
                                           treefile=trfn,
                                           schema_trf=schema_trf,
                                           otu_json=spInfoDict,
                                           ingroup_mrca=ingroup_mrca)

        # Prune sequnces below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()

        ids = IdDicts(conf, workdir=workdir)

        # Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape = FilterBlast(data_obj, ids)
        filteredScrape.blacklist = blacklist
        debug(add_unpubl_seq)
        if add_unpubl_seq is not None:
            filteredScrape.unpublished = True
            debug(filteredScrape.unpublished)
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.unpublished = True
            #filteredScrape.config.e_value_thresh = config['blast']['e_value_thresh']
            filteredScrape.backbone = True
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            filteredScrape.run_blast()
            debug(id_to_spn_addseq_json)
            filteredScrape.local_otu_json = id_to_spn_addseq_json
            filteredScrape.read_blast()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False
        else:
            # run the analysis
            sys.stdout.write("BLASTing input sequences\n")
            print(shared_blast_folder)
            if shared_blast_folder:
                filteredScrape.blast_subdir = shared_blast_folder
            else:
                shared_blast_folder = None
            filteredScrape.run_blast()
            filteredScrape.read_blast(blast_dir=shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            filteredScrape.dump()
            debug(threshold)
            sys.stdout.write("Filter the sequences\n")
            if threshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            sys.stdout.write("Calculate the phylogeny\n")
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
        filteredScrape.read_blast(blast_dir=shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        sys.stdout.write("Filter the sequences\n")
        if threshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.data.prune_short(0.75)
        sys.stdout.write("calculate the phylogeny\n")
        # with open("debug_PS_instances_w1", "a") as PS:
        #     PS.write("{}".format(filteredScrape.__dict__))
        # with open("debug_ATT_instances_w1", "a") as ATT:
        #     ATT.write("{}".format(filteredScrape.data.__dict__))
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()

    filteredScrape.write_otu_info(downtorank)
    filteredScrape.print_sp_d_recalc(downtorank)
    filteredScrape.print_sp_d_as_is()
    # with open("debug_PS_instances_f", "a") as PS:
    #     PS.write("{}".format(filteredScrape.__dict__))
    # with open("debug_ATT_instances_f", "a") as ATT:
    #     ATT.write("{}".format(filteredScrape.data.__dict__))
    return filteredScrape


def filter_data_run(seqaln,
                    mattype,
                    trfn,
                    schema_trf,
                    workdir,
                    threshold,
                    spInfoDict,
                    configfi,
                    selectby="blast",
                    downtorank=None,
                    blacklist=None,
                    add_unpubl_seq=None,
                    id_to_spn_addseq_json=None,
                    ingroup_mrca=None,
                    shared_blast_folder=None):
    """looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.
    """
    debug('Debugging mode is on')

    # debug(shared_blast_folder)
    # debug(some)
    # if _DEBUG_MK == 1:
    #     random.seed(3269235691)
    print(workdir)
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
        filteredScrape.repeat = 1   
    else:   
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        conf = ConfigObj(configfi)
        # print("config")
        debug(dir(conf))
        debug(conf.email)

        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                           mattype=mattype,
                                           workdir=workdir,
                                           treefile=trfn,
                                           schema_trf=schema_trf,
                                           otu_json=spInfoDict,
                                           ingroup_mrca=ingroup_mrca)

        # Prune sequnces below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()

        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()

        ids = IdDicts(conf, workdir=workdir)

        # Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape = FilterBlast(data_obj, ids)
        filteredScrape.blacklist = blacklist
        debug(add_unpubl_seq)
        if add_unpubl_seq is not None:
            filteredScrape.unpublished = True
            debug(filteredScrape.unpublished)
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.unpublished = True
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            filteredScrape.run_blast()
            debug(id_to_spn_addseq_json)
            filteredScrape.local_otu_json = id_to_spn_addseq_json
            filteredScrape.read_blast()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False
        else:
            # run the analysis
            sys.stdout.write("BLASTing input sequences\n")
            # uncomment next line if you want to have a shared blast folder and change the path to something meaningful. Remember to change the gifilename setting in the config file to true.
            print(shared_blast_folder)
            if shared_blast_folder:
                filteredScrape.blast_subdir = shared_blast_folder
            else:
                shared_blast_folder = None
            filteredScrape.run_blast()
            filteredScrape.read_blast(blast_dir=shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            filteredScrape.dump()
            debug(threshold)
            sys.stdout.write("Filter the sequences\n")
            if threshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            sys.stdout.write("Calculate the phylogeny\n")
            # with open("debug_PS_instances_r1", "a") as PS:
            #     PS.write("{}".format(filteredScrape.__dict__))
            # with open("debug_ATT_instances_r1", "a") as ATT:
            #     ATT.write("{}".format(filteredScrape.data.__dict__))
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
        filteredScrape.read_blast(blast_dir=shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        sys.stdout.write("Filter the sequences\n")
        if threshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.data.prune_short(0.75)
        sys.stdout.write("calculate the phylogeny\n")
        # with open("debug_PS_instances_w1", "a") as PS:
        #     PS.write("{}".format(filteredScrape.__dict__))
        # with open("debug_ATT_instances_w1", "a") as ATT:
        #     ATT.write("{}".format(filteredScrape.data.__dict__))
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()

    filteredScrape.write_otu_info(downtorank)
    filteredScrape.print_sp_d_recalc(downtorank)
    filteredScrape.print_sp_d_as_is()
    # with open("debug_PS_instances_f", "a") as PS:
    #     PS.write("{}".format(filteredScrape.__dict__))
    # with open("debug_ATT_instances_f", "a") as ATT:
    #     ATT.write("{}".format(filteredScrape.data.__dict__))
    return filteredScrape


# # # # # # # # # # # # # # # # # # # # # # #
def make_settings_class(seqaln, mattype, trfn, schema_trf, workdir, 
                        threshold=None, selectby=None, downtorank=None, spInfoDict=None, add_unpubl_seq=None, 
                        id_to_spn_addseq_json=None, configfi=None, blacklist=None, shared_blast_folder=None,
                        delay=None, trim=None):
    """all the settings are set here and can then be fed to the FilterClass
    """
    settings = Settings(seqaln=seqaln, mattype=mattype, trfn=trfn, schema_trf=schema_trf, workdir=workdir,
                        threshold=threshold, selectby=selectby, downtorank=downtorank, spInfoDict=spInfoDict,
                        add_unpubl_seq=add_unpubl_seq, id_to_spn_addseq_json=id_to_spn_addseq_json, configfi=configfi,
                        blacklist=blacklist, shared_blast_folder=shared_blast_folder, delay=delay, trim=trim)
    return settings
    

def run_with_settings(settings):
    """looks for pickeled file to continue run, or builds and runs
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output."""
    debug('Debugging mode is on')

    # if _DEBUG_MK == 1:
    #     random.seed(3269235691)
    if os.path.isfile("{}/scrape_checkpoint.p".format(settings.workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(settings.workdir), 'rb'))
        filteredScrape.repeat = 1   
    else: 
        conf = ConfigObj(settings.configfi)
        # print("config")
        debug(dir(conf))
        debug(conf.email)

        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=settings.seqaln, 
                                           mattype=settings.mattype,
                                           workdir=settings.workdir,
                                           treefile=settings.trfn,
                                           schema_trf=settings.schema_trf,
                                           otu_json=settings.spInfoDict,
                                           ingroup_mrca=None)

        # Prune sequences below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()

        data_obj.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()
        
        ids = IdDicts(conf, workdir=settings.workdir)

        filteredScrape = FilterBlast(data_obj, ids, settings)
        filteredScrape.write_otu_info(settings.downtorank)
        
        if settings.add_unpubl_seq is not None:
            filteredScrape.unpublished = True
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.write_unpubl_blastdb(settings.add_unpubl_seq)
            filteredScrape.run_blast(settings.delay)
            filteredScrape.local_otu_json = settings.id_to_spn_addseq_json
            filteredScrape.read_blast()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False

        # run the ananlyses
        if filteredScrape.unpublished is not True:
            filteredScrape.run_blast(settings.delay)
            filteredScrape.read_blast(blast_dir=settings.shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            filteredScrape.dump()
            if settings.threshold is not None:
                filteredScrape.sp_dict(settings.downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(threshold=settings.threshold, selectby=settings.selectby)
                filteredScrape.replace_new_seq()
            debug("from replace to streamed aln")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
    while filteredScrape.repeat is 1:
        filteredScrape.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
        filteredScrape.data.write_otus("otu_info", schema='table')
        filteredScrape.run_blast(settings.delay)
        filteredScrape.read_blast(blast_dir=settings.shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        if settings.threshold is not None:
            filteredScrape.sp_dict(settings.downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=settings.threshold, selectby=settings.selectby)
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

