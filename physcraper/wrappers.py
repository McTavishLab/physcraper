#!/usr/bin/env python
from __future__ import absolute_import
import pickle
import sys
import os
import shutil
import json
from physcraper import (
    generate_ATT_from_files,
    ConfigObj,
    IdDicts,
    PhyscraperScrape,
    AlignTreeTax,
    OtuJsonDict,
    debug
)
from physcraper.filterblast import FilterBlast
from dendropy import DnaCharacterMatrix

import physcraper.writeinfofiles as writeinfofiles
from physcraper.concat import Concat

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
#
#
# # TODO: not used, process_ott.sh does not exist
# def sync_ott(configfi):
#     conf = ConfigObj(configfi)
#     subprocess.call(["process_ott.sh", "".format(conf.ott_ncbi)])
#
#
# # TODO: not used
# # generates IdDicts physcrapper class
# def get_ottid(configfi, cwd):
#     conf = ConfigObj(configfi)
#     ids = IdDicts(conf, cwd)
#     return ids


def license_print():
    sys.stdout.write(
    """
    Physcraper: automatic updating of phylogenies
    Copyright (C) 2019  E.J. McTavish and M. Kandziora

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    """)


def load_ids_obj(conf, workdir):
    """
    Generates the IdDict class object.

    :param conf: Config Object of physcraper class
    :param workdir: working directory
    :return:
    """
    if os.path.isfile("{}/id_pickle.p".format(workdir)):
        sys.stdout.write("Reloading id dicts from {}\n".format(workdir))
        ids = pickle.load(open("{}/id_pickle.p".format(workdir), "rb"))
    else:
        sys.stdout.write("setting up ID dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, "{}/id_pickle.p".format(workdir))
        ids.dump("{}/id_pickle.p".format(workdir))
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
        # This is particularly important when using loci that have been de-concatenated,
        # as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label="^ot:ottTaxonName")
        data_obj.write_otus("otu_info", schema="table")
        data_obj.dump()
    assert isinstance(data_obj, AlignTreeTax)
    return data_obj


def make_otujsondict(id_to_spn, workdir, ids, local=False):
    """
    Generate a dictionary equivalent to the OToL one.

    :param id_to_spn: csv delimited file, where tipnames correspond to species names
    :param workdir: the working directory
    :param ids: physcraper Id object
    :param local: is needed for local database to change file name of otujson dict
    :return: otu dict as json file
    """
    workdir = os.path.abspath(workdir)
    # print(workdir)
    otu_jsonfi = "{}/otu_dict.json".format(workdir)
    if local is not False:
        otu_jsonfi = "{}/otu_dict_localseq.json".format(workdir)

    # print(id_to_spn)
    if os.path.exists(otu_jsonfi):
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = OtuJsonDict(id_to_spn, ids)
        # print(otu_json)
        # json.dump(otu_json, open("otu_dict.json", "w"))
        with open('{}'.format(otu_jsonfi), 'wb') as outfile:
            json.dump(otu_json, outfile)
        
        # json.dump(otu_json, open(otu_jsonfi, 'wb'))


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
        scraper = PhyscraperScrape(data_obj, ids, ingroup_mrca)
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
        write_out_files(scraper)
    writeinfofiles.get_additional_GB_info(scraper)
    return scraper


def PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids, selectby,
                  shared_blast_folder, threshold, ingroup_mrca, backbone=None):
    """
    This is the filtering mode for a Physcraper run:
    update aln and tre as long as new seqs are found, but filters the found sequences according to user settings.

    :param backbone: use existing tree as backbone during tree recalculation
    :param add_unpubl_seq: path to database if you want to add sequences from a local database
    :param blacklist: list with accession numbers of sequences not to be included
    :param data_obj: PS ATT object
    :param downtorank: delimits the rank for filtering, e.g. species/genus/ subfamily
    :param id_to_spn_addseq_json: JSON file where tip names correspond to species names for the local database
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
        filteredScrape = FilterBlast(data_obj, ids, ingroup_mrca)
        if backbone is True:
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
            filteredScrape.data.unpubl_otu_json = json.load(open("{}/otu_dict_localseq.json".format(data_obj.workdir)))
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            filteredScrape.run_blast_wrapper()
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
                if len(filteredScrape.new_seqs_otu_id) > 0:
                    filteredScrape.sp_dict(downtorank)
                    filteredScrape.make_sp_seq_dict()
                    filteredScrape.how_many_sp_to_keep(selectby=selectby)
                    filteredScrape.replace_new_seq()
            sys.stdout.write("Calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
            filteredScrape.data.write_otus("otu_info", schema="table")
            write_out_files(filteredScrape, downtorank)
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
            if len(filteredScrape.new_seqs_otu_id) > 0:
                filteredScrape.sp_dict(downtorank)
                filteredScrape.make_sp_seq_dict()
                filteredScrape.how_many_sp_to_keep(selectby=selectby)
                filteredScrape.replace_new_seq()
        filteredScrape.data.prune_short()
        sys.stdout.write("calculate the phylogeny\n")
        filteredScrape.generate_streamed_alignment()
        filteredScrape.dump()
        write_out_files(filteredScrape, downtorank)
        if backbone:
            filteredScrape.repeat = 0
    writeinfofiles.get_additional_GB_info(filteredScrape)
    return filteredScrape


def add_different_rank(seqaln,
                       mattype,
                       trfn,
                       schema_trf,
                       workdir,
                       threshold,
                       id_to_spn,
                       new_confifi,
                       selectby="blast",
                       downtorank=None,
                       blacklist=None,
                       add_unpubl_seq=None,
                       id_to_spn_addseq_json=None,
                       ingroup_mrca=None,
                       shared_blast_folder=None,
                       backbone=False):
    """looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found. 
    This uses the FilterBlast subclass to be able to filter the blast output.
    """
    license_print()
    debug("Debugging mode is on")

    dump_fn = "add_different_rank{}_{}.run".format(ingroup_mrca, downtorank)
    # if files does not exists, this loop was not yet run, if exitsts, go to next
    if os.path.isfile("{}/{}".format(workdir, dump_fn)):
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))
    else:

        assert os.path.isfile("{}/scrape_checkpoint.p".format(workdir))

        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir), 'rb'))

        # copy previous files to different folder
        count = 1
        while os.path.exists("{}/update_{}".format(workdir, count)):
            count += 1
        os.mkdir("{}/update_{}".format(workdir, count))
        old_runs = "{}/update_{}".format(workdir, count)

        src_files = os.listdir(workdir)
        for file_name in src_files:
            full_file_name = os.path.join(workdir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, old_runs)

        filteredScrape.repeat = 1
        conf = ConfigObj(new_confifi)
        # add new config
        assert filteredScrape.config != conf
        filteredScrape.config = conf
        assert filteredScrape.config == conf

        # set new ingroup_mrca
        filteredScrape.data.ott_mrca = ingroup_mrca
        filteredScrape.mrca_ncbi = filteredScrape.ids.ott_to_ncbi[filteredScrape.data.ott_mrca]
        assert filteredScrape.data.ott_mrca == ingroup_mrca

        with open(filteredScrape.logfile, "a") as log:
                log.write("You run 'add_different_rank' with the following settings: rank: {} and ingroup_mrca: {}. \n".format(downtorank, ingroup_mrca))

        # here the filter standard function continues...
        if backbone is True:
            filteredScrape.backbone = backbone
            filteredScrape.data.write_files(treepath="backbone.tre", alnpath="backbone.fas")
        else:
            filteredScrape.backbone = False
        # set new downtorank and numbers:
        filteredScrape.add_setting_to_self(downtorank, threshold)
        filteredScrape.blacklist = blacklist
    
        if add_unpubl_seq is not None:
            filteredScrape.unpublished = True
        if filteredScrape.unpublished is True:  # use unpublished data
            sys.stdout.write("Blasting against local unpublished data")
            filteredScrape.data.unpubl_otu_json = json.load(open("{}/otu_dict_localseq.json".format(workdir)))
            filteredScrape.write_unpubl_blastdb(add_unpubl_seq)
            filteredScrape.run_blast_wrapper()
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
                if len(filteredScrape.new_seqs_otu_id) > 0:
                    filteredScrape.sp_dict(downtorank)
                    filteredScrape.make_sp_seq_dict()
                    filteredScrape.how_many_sp_to_keep(selectby=selectby)
                    filteredScrape.replace_new_seq()
            sys.stdout.write("Calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
            filteredScrape.data.write_otus("otu_info", schema="table")
            write_out_files(filteredScrape, downtorank)
            if backbone:
                filteredScrape.repeat = 0
            # set back to normal - only used to reassess formerly discarded seq in first round
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
                if len(filteredScrape.new_seqs_otu_id) > 0:
                    filteredScrape.sp_dict(downtorank)
                    filteredScrape.make_sp_seq_dict()
                    filteredScrape.how_many_sp_to_keep(selectby=selectby)
                    filteredScrape.replace_new_seq()
            filteredScrape.data.prune_short()
            sys.stdout.write("calculate the phylogeny\n")
            filteredScrape.generate_streamed_alignment()
            filteredScrape.dump()
            write_out_files(filteredScrape, downtorank)
            if backbone:
                filteredScrape.repeat = 0
        writeinfofiles.get_additional_GB_info(filteredScrape)
        filteredScrape.dump()
    dump_fn = "add_different_rank{}_{}.run".format(ingroup_mrca, downtorank)
    fn = open(dump_fn, "w")
    fn.write("add different rank with following settings {} and {} finished".format(ingroup_mrca, downtorank))
    fn.close()
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
    license_print()
    debug("Debugging mode is on")
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    conf = ConfigObj(configfi)
    data_obj = load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir)
    # Mapping identifiers between OpenTree and NCBI requires an identifier dict object
    ids = load_ids_obj(conf, workdir, ingroup_mrca)
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
    license_print()
    debug("Debugging mode is on")
    if not os.path.exists(workdir):
        os.mkdir(workdir)
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
    license_print()

    debug("Debugging mode is on")
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    # read the config file into a configuration object
    conf = ConfigObj(configfi)
    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_otol_data(conf, ingroup_mrca, mattype, seqaln, study_id, tree_id, workdir)
    ids = load_ids_obj(conf, workdir)

    # make json file for unpublished database
    if add_unpubl_seq is not None:
        make_otujsondict(id_to_spn_addseq_json, workdir, ids, local=True)

    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold, ingroup_mrca)
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
    license_print()
    debug("Debugging mode is on")
    print(workdir)
    print(os.path.exists(workdir))
    if not os.path.exists(workdir):
        print("make wd")
        os.makedirs(workdir)
    conf = ConfigObj(configfi)
    ids = load_ids_obj(conf, workdir)

    make_otujsondict(id_to_spn, workdir, ids)
    # make json file for unpublished database
    if add_unpubl_seq is not None:
        make_otujsondict(id_to_spn_addseq_json, workdir, ids, local=True)

    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold, ingroup_mrca)
    save_copy_code(workdir)
    return filteredScrape


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
    license_print()

    # read the config file into a configuration object
    conf = ConfigObj(configfi)

    # Generate an linked Alignment-Tree-Taxa object
    data_obj = load_own_data(conf, seqaln, mattype, trfn, schema_trf, workdir, ingroup_mrca)
    ids = load_ids_obj(conf, workdir)
    filteredScrape = PS_filter_run(add_unpubl_seq, blacklist, data_obj, downtorank, id_to_spn_addseq_json, ids,
                                   selectby, shared_blast_folder, threshold, ingroup_mrca, backbone=True)
    save_copy_code(workdir)
    return filteredScrape


def concat(genelistdict, workdir_comb, email, num_threads=None, percentage=0.37, user_concat_fn=None, backbone=False):
    """This is to concatenate different physcraper runs into a single alignment and tree.
    genelistdict is a dict with gene names as key and the corresponding workdir
    """
    license_print()

    if not os.path.exists(path="{}/concat_checkpoint.p".format(workdir_comb)):
        if not os.path.exists(path="{}/load_single_data.p".format(workdir_comb)):
            # save_copy_code(workdir_comb)
            conc = Concat(workdir_comb, email)
            conc.concatfile = user_concat_fn
            for item in genelistdict.keys():
                conc.load_single_genes(genelistdict[item]["workdir"], genelistdict[item]["pickle"], item)
            conc.combine()
        else:
            sys.stdout.write("load single data dump file\n")
            conc = pickle.load(open("{}/load_single_data.p".format(workdir_comb), "rb"))
            # conc.dump()
        conc.sp_seq_counter()
        conc.get_largest_tre()
        conc.make_sp_gene_dict()
        conc.make_alns_dict()
        conc.concatenate_alns()
        conc.get_short_seq_from_concat(percentage)
        conc.remove_short_seq()
        conc.dump()
    else:
        sys.stdout.write("load concat_checkpoint dump file\n")
        conc = pickle.load(open("{}/concat_checkpoint.p".format(workdir_comb), "rb")) 
    conc.backbone = backbone
    conc.make_concat_table()
    conc.write_partition()
    conc.write_otu_info()
    conc.place_new_seqs(num_threads)
    
    if backbone is False:
        conc.calculate_bootstrap(num_threads)
        conc.write_labelled('RAxML_bestTree.autoMRE_fa')
    else:
        conc.est_full_tree(num_threads)
        conc.write_labelled('RAxML_bestTree.backbone_concat')
    return conc


def save_copy_code(workdir_comb):
    i = 1
    new_dir = "{}/physcraper_runcopy{}".format(workdir_comb, i)
    if os.path.exists("{}/physcraper_runcopy{}".format(workdir_comb, i)):
        while os.path.exists("{}/physcraper_runcopy{}".format(workdir_comb, i)):
            i += 1
        new_dir = "{}/physcraper_runcopy{}".format(workdir_comb, i)
    shutil.copytree("./physcraper", new_dir)


def write_out_files(obj, downtorank=None):
    """Wrapper function for writing information output files.

    Writes different output tables to file: Makes reading important information less code heavy.

    1. table with taxon names and sampling.
    2. a file with all relevant GenBank info to file (otu_dict).

    It uses the self.sp_d to get sampling information, that's why the downtorank is required.

    :param obj: either FilterBlast or PhyScraper Scrape object
    :param downtorank: hierarchical filter
    :return: writes output to file
    """

    writeinfofiles.write_otu_info(obj)
    if isinstance(obj, FilterBlast):
        writeinfofiles.taxon_sampling(obj, downtorank)


