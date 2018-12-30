#!/usr/bin/env python
from __future__ import absolute_import
import pickle
import sys
import os
import subprocess
import shutil
import json
from physcraper import (
    generate_ATT_from_phylesystem,
    generate_ATT_from_files,
    ConfigObj,
    IdDicts,
    PhyscraperScrape,
    AlignTreeTax,
    OtuJsonDict
)
from physcraper import FilterBlast, Settings, debug  # Concat
from dendropy import DnaCharacterMatrix
from .concat import Concat

print("Current Wrapper Version number: 12192018.0")

# TODO: we never do anything with the function nor the file
# def sync_ncbi(configfi):
#     conf = ConfigObj(configfi)
#     subprocess.call(
#         [
#             "rsync",
#             "av",
#             "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz",
#             "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp),
#         ]
#     )
#     subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


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


def load_ids_obj(conf, workdir):
    """
    Generates the IdDict class object.

    :param conf: Config Object of physcraper class
    :param workdir: working directory
    :return:
    """
    if os.path.isfile(conf.id_pickle):
        sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
        ids = pickle.load(open(conf.id_pickle, "rb"))
    else:
        sys.stdout.write("setting up ID dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        ids.dump(workdir)
    return ids


def load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir):
    """
    Generates ATT object from OToL data.

    :param conf: conf object from physcraper
    :param ingroup_mrca: mrca of ingroup as OTT ID
    :param mattype: alignment matrix type
    :param seqaln: alignment file name
    :param study_id: OToL study ID
    :param tree_id: OToL tree ID
    :param workdir: working directory
    :return: ATT object
    """
    if os.path.isfile("{}/att_checkpoint.p".format(workdir)):
        sys.stdout.write("Reloading data object from pickle file\n")
        data_obj = pickle.load(open("{}/att_checkpoint.p".format(workdir), "rb"))
    else:
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()

        aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_phylesystem(aln=aln,
                                                 workdir=workdir,
                                                 config_obj=conf,
                                                 study_id=study_id,
                                                 tree_id=tree_id,
                                                 phylesystem_loc=conf.phylesystem_loc,
                                                 ingroup_mrca=ingroup_mrca)
        # Prune sequences below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label="^ot:ottTaxonName")
        data_obj.write_otus("otu_info", schema="table")
        data_obj.dump()
    assert isinstance(data_obj, AlignTreeTax)
    return data_obj


def make_otujsondict(id_to_spn, workdir, ids):
    """
    Generate a dictionary equivalent to the OToL one.

    :param id_to_spn: csv delimited file, where tipnames correspond to species names
    :param ids: physcraper Id object
    :return: otu dict as json file
    """
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    if os.path.exists(otu_jsonfi):
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = OtuJsonDict(id_to_spn, ids)
        json.dump(otu_json, open(otu_jsonfi, "w"))


def load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca):
    """
    Generates ATT object from own data.

    :param conf: conf object from physcraper
    :param seqaln: sequence alignment file
    :param mattype: format of sequence alignment
    :param trfn: tree file
    :param schema_trf: format of tree file
    :param workdir: working directory
    :param ingroup_mrca: mrca of ingroup as OTT ID
    :return: ATT object
    """
    otu_jsonfi = "{}/otu_dict.json".format(workdir)
    assert os.path.exists(otu_jsonfi)

    if os.path.isfile("{}/att_checkpoint.p".format(workdir)):
        sys.stdout.write("Reloading from pickled scrapefile: ATT\n")
        data_obj = pickle.load(open("{}/att_checkpoint.p".format(workdir), "rb"))
    else:
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        # read the config file into a configuration object
        # Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln,
                                           mattype=mattype,
                                           workdir=workdir,
                                           config_obj=conf,
                                           treefile=trfn,
                                           schema_trf=schema_trf,
                                           otu_json=otu_jsonfi,
                                           ingroup_mrca=ingroup_mrca)

        # Prune sequences below a certain length threshold
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label="^ot:ottTaxonName")
        data_obj.write_otus("otu_info", schema="table")
        data_obj.dump()
    assert isinstance(data_obj, AlignTreeTax)
    return data_obj


def PS_standard_run(data_obj, ids, shared_blast_folder):
    """
    This is the standard mode for a Physcraper run:
    update aln and tre as long as new seqs are found, no filtering.

    :param data_obj: ATT object
    :param ids: IdDict object
    :param shared_blast_folder: path to folder for shared blast runs
    :return: PS run
    """
    if os.path.isfile("{}/scrape_checkpoint.p".format(data_obj.workdir)):
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        scraper = pickle.load(open("{}/scrape_checkpoint.p".format(data_obj.workdir), 'rb'))
        scraper.repeat = 1
    else:
        scraper = PhyscraperScrape(data_obj, ids)
        # run the analyses
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.run_blast_wrapper()
        scraper.read_blast_wrapper(blast_dir=shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
        scraper.dump("scrape_checkpoint.p")
    while scraper.repeat == 1:
        scraper.data.write_labelled(label="^ot:ottTaxonName")
        scraper.data.write_otus("otu_info", schema="table")
        if shared_blast_folder:
            scraper.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        scraper.run_blast_wrapper()
        scraper.read_blast_wrapper(blast_dir=shared_blast_folder)
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
        scraper.dump()
        scraper.write_otu_info()
    # scraper.write_otu_info()
    scraper.get_additional_GB_info()
    return scraper


def PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids, selectby,
                  shared_blast_folder, threshold, backbone=None):
    """
    This is the filtering mode for a Physcraper run:
    update aln and tre as long as new seqs are found, but filters the found sequences according to user settings.

    :param backbone: use existing tree as backbone during tree recalculation
    :param add_unpubl_seq: path to database if you want to add sequences from a local database
    :param blacklist: list with accession numbers of sequences not to be included
    :param data_obj: PS ATT object
    :param downtorank: delimits the rank for filtering, e.g. species/genus/ subfamily
    :param id_to_spn_addseq_json: JSON file whre tip names correspond to species names for the locl databse
    :param ids: PS IdDict object
    :param selectby: mode of filtering, either blast or length
    :param shared_blast_folder: path to folder for shared blast runs
    :param threshold: integer delimiting the number of sequences per OTU
    :return: PS run
    """
    if os.path.isfile("{}/scrape_checkpoint.p".format(data_obj.workdir)):
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(data_obj.workdir), 'rb'))
        filteredScrape.repeat = 1
    else:
        # Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape = FilterBlast(data_obj, ids)
        if backbone == True:
            filteredScrape.backbone = backbone
            filteredScrape.data.write_files(treepath="backbone.tre", alnpath="backbone.fas")
        else:
            filteredScrape.backbone = False
        filteredScrape.add_setting_to_self(downtorank, threshold)
        filteredScrape.blacklist = blacklist
        if add_unpubl_seq is not None:
            filteredScrape.unpublished = True
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.unpublished = True
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            filteredScrape.run_blast_wrapper()
            filteredScrape.data.unpubl_otu_json = id_to_spn_addseq_json
            filteredScrape.read_blast_wrapper()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False
            if backbone:
                filteredScrape.repeat = 1
        else:
            sys.stdout.write("BLASTing input sequences\n")
            if shared_blast_folder:
                filteredScrape.blast_subdir = shared_blast_folder
            else:
                shared_blast_folder = None
            filteredScrape.run_blast_wrapper()
            filteredScrape.read_blast_wrapper(blast_dir=shared_blast_folder)
            filteredScrape.remove_identical_seqs()
            sys.stdout.write("Filter the sequences\n")
            if threshold is not None:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
                filteredScrape.replace_new_seq()
            sys.stdout.write("Calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
            filteredScrape.data.write_otus("otu_info", schema="table")
            filteredScrape.write_out_files(downtorank)
            if backbone:
                filteredScrape.repeat = 0
    while filteredScrape.repeat == 1:
        filteredScrape.data.write_labelled(label="^ot:ottTaxonName", add_gb_id=True)
        filteredScrape.data.write_otus("otu_info", schema="table")
        sys.stdout.write("BLASTing input sequences\n")
        if shared_blast_folder:
            filteredScrape.blast_subdir = shared_blast_folder
        else:
            shared_blast_folder = None
        filteredScrape.run_blast_wrapper()
        filteredScrape.read_blast_wrapper(blast_dir=shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        sys.stdout.write("Filter the sequences\n")
        if threshold is not None:
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=threshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.data.prune_short()
        sys.stdout.write("calculate the phylogeny\n")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        filteredScrape.write_otu_info()
        filteredScrape.write_out_files(downtorank)
        if backbone:
            filteredScrape.repeat = 0
    filteredScrape.write_out_files(downtorank)
    filteredScrape.get_additional_GB_info()
    return filteredScrape

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
    debug("Debugging mode is on")
    conf = ConfigObj(configfi)
    data_obj = load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir)
    # Mapping identifiers between OpenTree and NCBI requires an identifier dict object
    ids = load_ids_obj(conf, workdir)
    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    scraper = PS_standard_run(data_obj, ids, shared_blast_folder)
    save_copy_code(workdir)
    return scraper


def own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 id_to_spn,
                 configfi,
                 ingroup_mrca=None,
                 shared_blast_folder=None):
    """This is the wrapper function to start a PhyScraper standard run with your own data.
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

    debug("Debugging mode is on")
    conf = ConfigObj(configfi)
    ids = load_ids_obj(conf, workdir)

    make_otujsondict(id_to_spn, workdir, ids)
    data_obj = load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
    # Mapping identifiers between original data and NCBI requires an identifier dict object
    # scraper = PhyscraperScrape(data_obj, ids)
    scraper = PS_standard_run(data_obj, ids, shared_blast_folder)
    save_copy_code(workdir)
    return 1


def filter_OTOL(study_id,
                tree_id,
                seqaln,
                mattype,
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

    This uses the FilterBlast subclass to be able to filter the blast output using data from OToL."""
    debug("Debugging mode is on")
    # read the config file into a configuration object
    conf = ConfigObj(configfi)
    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir)
    ids = load_ids_obj(conf, workdir)
    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold)
    save_copy_code(workdir)
    return filteredScrape


def filter_data_run(seqaln,
                    mattype,
                    trfn,
                    schema_trf,
                    workdir,
                    threshold,
                    id_to_spn,
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
    debug("Debugging mode is on")
    conf = ConfigObj(configfi)
    ids = load_ids_obj(conf, workdir)

    make_otujsondict(id_to_spn, workdir, ids)

    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold)
    save_copy_code(workdir)
    return filteredScrape

# # # # # # # # # # # # # #

def add_unpubl_to_backbone(seqaln,
                           mattype,
                           trfn,
                           schema_trf,
                           workdir,
                           sp_info_jsonfi,
                           configfi,
                           add_unpubl_seq,
                           id_to_spn_addseq_json,
                           selectby=None,
                           downtorank=None,
                           threshold=None,
                           blacklist=None,
                           ingroup_mrca=None,
                           shared_blast_folder=None):
    """
    This uses the FilterBlast subclass to be able to filter the blast output.
    It adds unpublished data to an input tree (evalue should be higher than usual).
    Backbone will not be updated
    """

    # read the config file into a configuration object
    conf = ConfigObj(configfi)

    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
    ids = load_ids_obj(conf, workdir)
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold, backbone=True)
    save_copy_code(workdir)
    return filteredScrape


def concat(genelistdict, workdir_comb, email, num_threads=None, percentage=0.37, user_concat_fn=None, backbone=None):
    """This is to concatenate different physcraper runs into a single alignment and tree.
    genelistdict is a dict with gene names as key and the corresponding workdir
    """
    concat = Concat(workdir_comb, email)
    concat.concatfile = user_concat_fn
    for item in genelistdict.keys():
        concat.load_single_genes(genelistdict[item]["workdir"], genelistdict[item]["pickle"], item)
    concat.combine()
    concat.sp_seq_counter()
    concat.get_largest_tre()
    concat.make_sp_gene_dict()
    concat.make_alns_dict()
    concat.concatenate_alns()
    concat.get_short_seq_from_concat(percentage)
    concat.remove_short_seq()
    concat.make_concat_table()
    concat.write_partition()
    concat.place_new_seqs(num_threads)
    concat.est_full_tree(num_threads)
    concat.calculate_bootstrap(num_threads)
    concat.write_otu_info()
    save_copy_code(workdir_comb)

    return concat


def save_copy_code(workdir_comb):
    i = 1
    if os.path.exists("{}/physcraper_runcopy".format(workdir_comb)):
        prev_dir = "{}/physcraper_runcopy{}".format(workdir_comb, i)
        i += 1
        shutil.move("{}/physcraper_runcopy".format(workdir_comb), prev_dir)
    shutil.copytree("./physcraper", "{}/physcraper_runcopy".format(workdir_comb))


# # # # # # # # # # # # # # # # # # # # # # #
def make_settings_class(seqaln, mattype, trfn, schema_trf, workdir,
                        threshold=None, selectby=None, downtorank=None, spInfoDict=None, add_unpubl_seq=None,
                        id_to_spn_addseq_json=None, configfi=None, blacklist=None, shared_blast_folder=None):
    """all the settings are set here and can then be fed to the FilterClass
    """
    settings = Settings(seqaln=seqaln, mattype=mattype, trfn=trfn, schema_trf=schema_trf, workdir=workdir,
                        threshold=threshold, selectby=selectby, downtorank=downtorank, spInfoDict=spInfoDict,
                        add_unpubl_seq=add_unpubl_seq, id_to_spn_addseq_json=id_to_spn_addseq_json, configfi=configfi,
                        blacklist=blacklist, shared_blast_folder=shared_blast_folder)
    return settings


def run_with_settings(settings):
    """looks for pickeled file to continue run, or builds and runs
    new analysis for as long as new seqs are found.
    This uses the FilterBlast subclass to be able to filter the blast output."""
    debug("Debugging mode is on")
    if os.path.isfile("{}/scrape_checkpoint.p".format(settings.workdir)):
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(
            open("{}/scrape_checkpoint.p".format(settings.workdir), "rb")
        )
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
                                           config_obj=conf,
                                           treefile=settings.trfn,
                                           schema_trf=settings.schema_trf,
                                           otu_json=settings.spInfoDict,
                                           ingroup_mrca=None)

        # Prune sequences below a certain length threshold
        # This is particularly important when using loci that have been de-concatenated,
        # as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()

        data_obj.write_labelled(label="^ot:ottTaxonName", add_gb_id=True)
        data_obj.write_otus("otu_info", schema="table")
        data_obj.dump()

        ids = IdDicts(conf, workdir=settings.workdir)

        filteredScrape = FilterBlast(data_obj, ids, settings)
        filteredScrape.add_setting_to_self(settings.downtorank, settings.threshold)

        filteredScrape.write_out_files(settings.downtorank)

        if settings.add_unpubl_seq is not None:
            filteredScrape.unpublished = True
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.write_unpubl_blastdb(settings.add_unpubl_seq)
            filteredScrape.run_blast_wrapper()
            filteredScrape.local_otu_json = settings.id_to_spn_addseq_json
            filteredScrape.read_blast_wrapper()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment()
            filteredScrape.unpublished = False

        # run the ananlyses
        if filteredScrape.unpublished is not True:
            filteredScrape.run_blast_wrapper()
            filteredScrape.read_blast_wrapper(blast_dir=settings.shared_blast_folder)
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
        filteredScrape.data.write_labelled(label="^ot:ottTaxonName", add_gb_id=True)
        filteredScrape.data.write_otus("otu_info", schema="table")
        filteredScrape.run_blast_wrapper()
        filteredScrape.read_blast_wrapper(blast_dir=settings.shared_blast_folder)
        filteredScrape.remove_identical_seqs()
        if settings.threshold is not None:
            filteredScrape.sp_dict(settings.downtorank)
            filteredScrape.make_sp_seq_dict()
            filteredScrape.how_many_sp_to_keep(threshold=settings.threshold, selectby=settings.selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        filteredScrape.write_out_files(settings.downtorank)
    filteredScrape.get_additional_GB_info()
    shutil.copytree("./physcraper", "{}/physcraper_runcopy".format(settings.workdir))

    return filteredScrape