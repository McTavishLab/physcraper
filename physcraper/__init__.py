#!/usr/bin/env python
"""Physcraper module"""

from __future__ import absolute_import

import sys
import re
import os
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import contextlib
import time
import csv
# from mpi4py import MPI
from past.builtins import xrange
from builtins import input
from copy import deepcopy
from ete2 import NCBITaxa
import physcraper.AWSWWW as AWSWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI, APIWrapper
from peyotl.sugar import tree_of_life, taxomachine
from peyotl.nexson_syntax import (
    extract_tree,
    get_subtree_otus,
    extract_otu_nexson,
    PhyloSchema
)

# extension functions
from . import concat  # is the local concat class
from . import ncbi_data_parser  # is the ncbi data parser class and associated functions
from . import filter_by_local_blast  # functions for the FilterBlast filtering
from . import opentree_helpers
from . import writeinfofiles

if sys.version_info < (3,):
    from urllib2 import HTTPError
else:
    from urllib.error import HTTPError

_DEBUG = 1
_DEBUG_MK = 1
_deep_debug = 0

_VERBOSE = 0


@contextlib.contextmanager
def cd(path):
    # print 'initially inside {0}'.format(os.getcwd())
    CWD = os.getcwd()
    os.chdir(path)
    # print 'inside {0}'.format(os.getcwd())
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        # print 'finally inside {0}'.format(os.getcwd())
        os.chdir(CWD)


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


def deep_debug(msg):
    """short debugging command
    """
    if _deep_debug == 1:
        print(msg)


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


# which python physcraper file do I use?
debug("Current --init-- version number: 12-17-2018.0")
debug(os.path.realpath(__file__))


def get_user_input():
    """Asks for yes or no user input.

    :return: user input
    """
    debug("get user input")
    is_valid = 0
    x = None
    while not is_valid:
        try:
            x = input("Please write either 'yes' or 'no': ")
            if x == "yes" or x == "no":
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError as e:
            print("'{}' is not a valid answer.".format(e.args[0].split(": ")[1]))
    return x


class ConfigObj(object):
    """
    To build the class the following is needed:

      * **configfi**: a configuration file in a specific format, e.g. to read in self.e_value_thresh.

        The file needs to have a heading of the format: [blast] and then somewhere below that heading a string e_value_thresh = value

      * **interactive**: defaults to True, is used to interactively update the local blast databases

    During the initializing process the following self objects are generated:

      * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
      * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
      * **self.seq_len_perc**: value from 0 to 1. Defines how much shorter new seq can be compared to input
      * **self.trim_perc**: value that determines how many seq need to be present before the beginning and end of alignment will be trimmed
      * **self.maxlen**: max length for values to add to aln
      * **self.get_ncbi_taxonomy**: Path to sh file doing something...
      * **self.phylesystem_loc**: defines which phylesystem for OpenTree datastore is used - default: api, but can run on local version too.
      * **self.ott_ncbi**: file containing OTT id, ncbi and taxon name (??)
      * **self.id_pickle**: path to pickle file
      * **self.email**: email address used for blast queries
      * **self.blast_loc**: defines which blasting method to use:

          * either web-query (=remote)
          * from a local blast database (=local)
      * **self.num_threads**: number of cores to be used during a run
      * **self.url_base**:

          * if blastloc == remote: it defines the url for the blast queries.
          * if blastloc == local: url_base = None
      * **self.unmapped**: used for OToL original tips that can not be assigned to a taxon

          * keep: keep the unmapped taxa and asign them to life
          * remove: remove the unmapped taxa from aln and tre
      * **self.delay**: defines when to reblast sequences in days
      * **self.add_lower_taxa**: T/F, enables to re-access formerly filtered seq by allowing them be passed into remove_identical. Used if we first filter for higher rank and then want to filter for a lower rank.
      * **optional self.objects**:

          * if blastloc == local:

              * self.blastdb: this defines the path to the local blast database
              * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
              * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
    """

    def __init__(self, configfi, interactive=None):
        debug(configfi)
        debug(os.path.isfile(configfi))
        if _DEBUG:
            sys.stdout.write("Building config object\n")
        assert os.path.isfile(configfi), "file `{}` does not exists".format(configfi)
        config = configparser.ConfigParser()
        config.read_file(open(configfi))

        # read in blast settings
        self.email = config["blast"]["Entrez.email"]
        assert "@" in self.email, "your email `{}` does not have an @ sign".format(self.email)

        self.e_value_thresh = config["blast"]["e_value_thresh"]
        assert is_number(self.e_value_thresh), (
                "value `{}` does not exists".format(self.e_value_thresh)
        )
        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), (
                "value `{}`is not a number".format(self.e_value_thresh)
        )
        self.blast_loc = config["blast"]["location"]
        assert self.blast_loc in ["local", "remote"], (
                "your blast location `{}` is not remote or local".format(self.email)
        )
        if self.blast_loc == "local":
            self.blastdb = config["blast"]["localblastdb"]
            self.url_base = None
            self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
            self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]
        if self.blast_loc == "remote":
            self.url_base = config["blast"].get("url_base")
        if _DEBUG:
            sys.stdout.write("{}\n".format(self.email))
            if self.blast_loc == "remote":
                sys.stdout.write("url base = {}\n".format(self.url_base))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == "local":
                sys.stdout.write("local blast db {}\n".format(self.blastdb))
        self.num_threads = config["blast"].get("num_threads")
        debug("slurm threads")
        debug(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))

        debug(self.num_threads)
        self.gb_id_filename = config["blast"].get("gb_id_filename", False)
        if self.gb_id_filename is not False:
            if self.gb_id_filename == "True" or self.gb_id_filename == "true":
                self.gb_id_filename = True
            else:
                self.gb_id_filename = False
        debug("shared blast folder? {}".format(self.gb_id_filename))
        self.delay = int(config["blast"]["delay"])
        assert is_number(self.delay), (
                "value `{}`is not a number".format(self.delay)
        )
        # #############
        # read in physcraper settings
        self.unmapped = config["physcraper"]["unmapped"]
        assert self.unmapped in ["remove", "keep"], (
                "your unmapped statement '{}'' in the config file is not remove or keep"
               .format(self.unmapped)
        )
        self.seq_len_perc = float(config["physcraper"]["seq_len_perc"])
        assert 0 < self.seq_len_perc <= 1, (
                "value `{}` is not between 0 and 1".format(self.seq_len_perc)
        )
        self.trim_perc = float(config["physcraper"]["trim_perc"])
        assert 0 < self.trim_perc < 1, (
                "value `{}` is not between 0 and 1".format(self.trim_perc)
        )
        self.maxlen = float(config["physcraper"]["max_len"])
        assert 1 < self.maxlen, (
                "value `{}` is not larger than 1".format(self.maxlen)
        )
        self.add_lower_taxa = config["physcraper"]["add_lower_taxa"]
        if self.add_lower_taxa == "True" or self.add_lower_taxa == "true":
            self.add_lower_taxa = True
        else:
            self.add_lower_taxa = False
        assert self.add_lower_taxa in [True, False], (
                "self.add_lower_taxa `{}` is not True or False".format(self.add_lower_taxa)
        )
        # read in settings for internal Physcraper processes
        # default is api, but can run on local version of OpenTree datastore
        self.phylesystem_loc = config["phylesystem"]["location"]
        assert self.phylesystem_loc in ["local", "api"], \
            (
                "phylesystem location must be either local or api")
        self.ott_ncbi = config["taxonomy"][
            "ott_ncbi"
        ]
        assert os.path.isfile(self.ott_ncbi), (
                "file `{}` does not exists".format(self.ott_ncbi)
        )
        # rewrites relative path to absolute path so that it behaves when changing dirs
        self.id_pickle = os.path.abspath(config["taxonomy"]["id_pickle"])

        ####
        # check database status
        if interactive is None:
            interactive = sys.stdin.isatty()
            if interactive is False:
                sys.stdout.write("REMEMBER TO UPDATE THE NCBI DATABASES REGULARLY!!\n")
        if interactive is True:
            self._download_ncbi_parser()
            self._download_localblastdb()
        debug("check db file status?: {}".format(interactive))

    def _download_localblastdb(self):
        """Check if files are present and if they are uptodate.
        If not files will be downloaded.
        """
        if self.blast_loc == "local":
            # next line of codes exists to have interactive mode enabled while testing
            # this allows to not actually have a local ncbi database downloaded
            if not os.path.isfile("{}/empty_local_db_for_testing.nhr".format(self.blastdb)):
                if not os.path.isfile("{}/nt.69.nhr".format(self.blastdb)):
                    print("Do you want to download the blast nt databases from ncbi? Note: "
                          "This is a US government website! You agree to their terms")
                    x = get_user_input()
                    if x == "yes":
                        os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'"
                                  "{}/".format(self.blastdb))
                        os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'"
                                  "{}/".format(self.blastdb))
                        with cd(self.blastdb):
                            os.system("update_blastdb nt")
                            os.system("cat *.tar.gz | tar -xvzf - -i")
                            os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                            os.system("rm *.tar.gz*")
                    elif x == "no":
                        print(
                            "You did not agree to download data from ncbi. Program will default to blast web-queries.")
                        self.blast_loc = "remote"
                    else:
                        print("You did not type yes or no!")
                else:
                    download_date = os.path.getmtime("{}/nt.60.nhr".format(self.blastdb))
                    download_date = datetime.datetime.fromtimestamp(download_date)
                    today = datetime.datetime.now()
                    time_passed = (today - download_date).days
                    if time_passed >= 90:
                        print("""Your databases might not be uptodate anymore. 
                              You downloaded them {} days ago. Do you want to update the blast databases from ncbi? 
                              Note: This is a US government website! You agree to their terms.""".format(time_passed))
                        x = get_user_input()
                        if x == "yes":
                            with cd(self.blastdb):
                                os.system("update_blastdb nt")
                                os.system("cat *.tar.gz | tar -xvzf - -i")
                                os.system("update_blastdb taxdb")
                                os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                                os.system("rm *.tar.gz*")
                        elif x == "no":
                            print("You did not agree to update data from ncbi. Old database files will be used.")
                        else:
                            print("You did not type 'yes' or 'no'!")

    def _download_ncbi_parser(self):
        """Check if files are present and if they are up to date.
        If not files will be downloaded.
        """
        if self.blast_loc == "local":
            if not os.path.isfile(self.ncbi_parser_nodes_fn):
                print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms")
                x = get_user_input()
                if x == "yes":
                    os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                    os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                    os.system("mv nodes.dmp ./tests/data/")
                    os.system("mv names.dmp ./tests/data/")
                    os.system("rm taxdump.tar.gz")
                elif x == "no":
                    print("You did not agree to download data from ncbi. Program will default to blast web-queries.")
                    print("This is slow and crashes regularly!")
                    self.blast_loc = "remote"
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime(self.ncbi_parser_nodes_fn)
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days
                if time_passed >= 90:
                    print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms")
                    x = get_user_input()
                    if x == "yes":
                        os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                        os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                        os.system("mv nodes.dmp ./tests/data/")
                        os.system("mv names.dmp ./tests/data/")
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type yes or no!")


def get_dataset_from_treebase(study_id, phylesystem_loc="api"):
    """Function is used to get the aln from treebase, for a tree that OpenTree has the mapped tree.
    """
    try:
        nexson = get_nexson(study_id, phylesystem_loc)
    except HTTPError as err:
        sys.stderr.write(err)
        sys.stderr.write("couldn't find study id {} in phylesystem location {}\n".format(study_id, phylesystem_loc))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treebase record associated with study ")
        sys.exit(-2)
    else:
        tb_id = treebase_url.split(':S')[1]
        url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexml".format(tb_id)
        if _DEBUG:
            sys.stderr.write(url + "\n")
        dna = DataSet.get(url=url, schema="nexml")
        return dna


# ATT is a dumb acronym for Alignment Tree Taxa object
def generate_ATT_from_phylesystem(aln,
                                  workdir,
                                  config_obj,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api',
                                  ingroup_mrca=None):
    """gathers together tree, alignment, and study info - forces names to otu_ids.

    Study and tree ID's can be obtained by using python ./scripts/find_trees.py LINEAGE_NAME

    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.

    :param aln: dendropy :class:`DnaCharacterMatrix <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix>` alignment object
    :param workdir: path to working directory
    :param config_obj: config class containing the settings
    :param study_id: OToL study id of the corresponding phylogeny which shall be updated
    :param tree_id: OToL corresponding tree ID as some studies have several phylogenies
    :param phylesystem_loc: access the github version of the OpenTree data store, or a local clone
    :param ingroup_mrca: optional.  OToL identifier of the mrca of the clade that shall be updated (can be subset of the phylogeny)
    :return: object of class ATT
    """
    assert isinstance(aln, datamodel.charmatrixmodel.DnaCharacterMatrix), \
            "your alignment `{}` ist not of type DnaCharacterMatrix".format(aln)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore
    nexson = get_nexson(study_id, phylesystem_loc)
    newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
    newick = newick.replace(" ", "_")  # UGH Very heavy handed, need to make sure happens on alignment side as well.
    tre = Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    # this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    otus = get_subtree_otus(nexson, tree_id=tree_id)
    otu_dict = {}
    orig_lab_to_otu = {}
    treed_taxa = {}
    for otu_id in otus:
        otu_dict[otu_id] = extract_otu_nexson(nexson, otu_id)[otu_id]
        otu_dict[otu_id]["^physcraper:status"] = "original"
        otu_dict[otu_id]["^physcraper:last_blasted"] = "1800/01/01"
        orig = otu_dict[otu_id].get(u"^ot:originalLabel").replace(" ", "_")
        orig_lab_to_otu[orig] = otu_id
        treed_taxa[orig] = otu_dict[otu_id].get(u"^ot:ottId")
    for tax in aln.taxon_namespace:
        if tax .label in otu_dict:
            sys.stdout.write("{} aligned\n".format(tax.label))
        else:
            try:
                tax.label = orig_lab_to_otu[tax.label].encode("ascii")
            except KeyError:
                sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignment. "
                                 "This may indicate a mismatch between tree and alignment\n".format(tax.label))
    # need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    ott_ids = get_subtree_otus(nexson,
                               tree_id=tree_id,
                               subtree_id="ingroup",
                               return_format="ottid")
    if ingroup_mrca:
        if type(ingroup_mrca) == list:
            ott_ids = set(ingroup_mrca)
            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(ingroup_mrca)
    elif ott_ids:  # if no ingroup is specified, ott_ids will be none
        ott_mrca = get_mrca_ott(ott_ids)
    else:  # just get the mrca for teh whole tree
        ott_mrca = get_mrca_ott([otu_dict[otu_id].get(u"^ot:ottId") for otu_id in otu_dict])
    workdir = os.path.abspath(workdir)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, config_obj=config_obj)
    # newick should be bare, but alignment should be DNACharacterMatrix


def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            config_obj,
                            treefile,
                            otu_json,
                            schema_trf,
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem, use your own files instead.

    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.

    Note: has test -> test_owndata.py

    :param seqaln: path to sequence alignment
    :param mattype: string containing format of sequence alignment
    :param workdir: path to working directory
    :param config_obj: config class including the settings
    :param treefile: path to phylogeny
    :param otu_json: path to json file containing the translation of tip names to taxon names, generated with OtuJsonDict()
    :param schema_trf: string defining the format of the input phylogeny
    :param ingroup_mrca: optional - OToL ID of the mrca of the clade of interest. If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.

    :return: object of class ATT
    """

    # replace ? in seqaln with - : papara handles them as different characters
    with open(seqaln, "r") as fin:
        filedata = fin.read()
    filedata = filedata.replace("?", "-")
    # Write the file out again
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    new_seq_file = "{}/replaced_inputaln.fasta".format(workdir)
    with open("{}/replaced_inputaln.fasta".format(workdir), "w") as aln_file:
        aln_file.write(filedata)
    # use replaced aln as input
    aln = DnaCharacterMatrix.get(path=new_seq_file, schema=mattype)
    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore
    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    assert tre.taxon_namespace is aln.taxon_namespace, "tre and aln have not the same namespace."
    otu_newick = tre.as_string(schema=schema_trf)
    otu_dict = json.load(open(otu_json, "r"))
   # debug("get mrca: ")
   # debug(ingroup_mrca)
    if ingroup_mrca:
        if type(ingroup_mrca) == list:
            ott_ids = set(ingroup_mrca)
            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(ingroup_mrca)
    else:
        ott_ids = [otu_dict[otu].get(u'^ot:ottId', ) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir,
                        config_obj=config_obj, schema=schema_trf)


def standardize_label(item):
    """Make sure that the tip names are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

    :param item: original tip name
    :return: tip name in unicode
    """
    item_edit = item.replace("-", "")
    item_edit = item_edit.replace(" ", "")
    item_edit = item_edit.replace("_", "")
    item_edit = item_edit.replace("'", "")
    item_edit = item_edit.replace("/", "")
    return item_edit


def get_ott_taxon_info(spp_name):
    """get ottid, taxon name, and ncbid (if present) from Open Tree Taxonomy.
    ONLY works with version 3 of Open tree APIs

    :param spp_name: species name
    :return:
    """
    #This is only used to write out the opentree info file. Could use NCBI id's instead of name, and likely be quicker.
    # debug(spp_name)
    try:
        res = taxomachine.TNRS(spp_name)["results"][0]
    except IndexError:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0
    if res['matches'][0]['is_approximate_match'] == 1:
        sys.stderr.write("""exact match to taxon {} not found in open tree taxonomy.
                          Check spelling. Maybe {}?""".format(spp_name, res['matches'][0][u'ot:ottTaxonName']))
        return 0
    if res["matches"][0]["is_approximate_match"] == 0:
        ottid = res["matches"][0]["taxon"][u"ott_id"]
        ottname = res["matches"][0]["taxon"][u"unique_name"]
        ncbi_id = None
        for source in res["matches"][0]["taxon"][u"tax_sources"]:
            if source.startswith("ncbi"):
                ncbi_id = source.split(":")[1]
        return ottid, ottname, ncbi_id
    else:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0


def OtuJsonDict(id_to_spn, id_dict):
    """Makes otu json dict, which is also produced within the openTreeLife-query.

     This function is used, if files that shall be updated are not part of the OpenTreeofLife project.
    It reads in the file that contains the tip names and the corresponding species names.
    It then tries to get the different identifier from the OToL project or if not from ncbi.

    Reads input file into the var sp_info_dict, translates using an IdDict object
    using web to call Open tree, then ncbi if not found.

    :param id_to_spn: user file, that contains tip name and corresponding sp name for input files.
    :param id_dict: uses the id_dict generated earlier
    :return: dictionary with key: "otu_tiplabel" and value is another dict with the keys '^ncbi:taxon',
                                                    '^ot:ottTaxonName', '^ot:ottId', '^ot:originalLabel',
                                                    '^user:TaxonName', '^physcraper:status', '^physcraper:last_blasted'
    """
    sys.stdout.write("Set up OtuJsonDict \n")
    sp_info_dict = {}
    nosp = []
    with open(id_to_spn, mode="r") as infile:
        for lin in infile:
            ottid, ottname, ncbiid = None, None, None
            tipname, species = lin.strip().split(",")
            clean_lab = standardize_label(tipname)
            assert clean_lab not in sp_info_dict, ("standardized label ('{}') of `{}` already exists".format(clean_lab, tipname))
            otu_id = "otu{}".format(clean_lab)
            spn = species.replace("_", " ")
            info = get_ott_taxon_info(spn)
            if info:
                ottid, ottname, ncbiid = info
            if not info:
                ncbi = NCBITaxa()
                name2taxid = ncbi.get_name_translator([spn])
                if len(name2taxid.items()) >= 1:
                    ncbiid = name2taxid.items()[0][1][0]
                else:
                    sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                     "Proceeding without taxon info\n".format(spn))
                    nosp.append(spn)
            ncbi_spn = None
            if ncbiid is not None:
                ncbi_spn = spn
            else:
                ncbi_spn = id_dict.ott_to_ncbi[ottid]
            sp_info_dict[otu_id] = {
                "^ncbi:taxon": ncbiid,

                "^ot:ottTaxonName": ottname,
                "^ot:ottId": ottid,
                "^ot:originalLabel": tipname,
                "^user:TaxonName": species,
                "^physcraper:status": "original",
                "^physcraper:last_blasted": "1900/01/01",
                }
            if ncbi_spn is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ncbi_spn
                sp_info_dict[otu_id]["^ncbi:TaxonName"] = ncbi_spn
            elif ottname is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ottname
            elif sp_info_dict[otu_id]['^user:TaxonName']:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = sp_info_dict[otu_id]['^user:TaxonName']
            assert sp_info_dict[otu_id]["^physcraper:TaxonName"]  # is not None
    return sp_info_dict


class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
        Hypothetically, all the keys in the  otu_dict should be clean.

        To build the class the following is needed:

          * **newick**: dendropy.tre.as_string(schema=schema_trf) object
          * **otu_dict**: json file including the otu_dict information generated earlier
          * **alignment**: dendropy :class:`DnaCharacterMatrix <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix>` object
          * **ingroup_mrca**: OToL identifier of the group of interest, either subclade as defined by user or of all tip labels in the phylogeny
          * **workdir**: the path to the corresponding working directory
          * **config_obj**: Config class
          * **schema**: optional argument to define tre file schema, if different from "newick"

        During the initializing process the following self objects are generated:

          * **self.aln**: contains the alignment and which will be updated during the run
          * **self.tre**: contains the phylogeny, which will be updated during the run
          * **self.otu_dict**: dictionary with taxon information and physcraper relevant stuff

               * key: a unique identifier (otu plus either "tiplabel of phylogeny" or for newly found sequences PS_number.
               * value: dictionary with the following key:values:

                    * '^ncbi:gi': GenBank identifier - deprecated by Genbank - only older sequences will have it
                    * '^ncbi:accession': Genbanks accession number
                    * '^ncbi:title': title of Genbank sequence submission
                    * '^ncbi:taxon': ncbi taxon identifier
                    * '^ot:ottId': OToL taxon identifier
                    * '^physcraper:status': contains information if it was 'original', 'queried', 'removed', 'added during filtering process'
                    * '^ot:ottTaxonName': OToL taxon name
                    * '^physcraper:last_blasted': contains the date when the sequence was blasted.

                         If the year is different from the 20th century, it tells us
                         something about the initial status:
                         * 1800 = never blasted, not yet considered to be added
                         * 1900 = never blasted and not added - see status for more information
                         * this century = blasted and added.
                    * '^user:TaxonName': optional, user given label from OtuJsonDict
                    * "^ot:originalLabel" optional, user given tip label of phylogeny
          * **self.ps_otu**: iterator for new otu IDs, is used as key for self.otu_dict
          * **self.workdir**: contains the path to the working directory, if folder does not exists it is generated.
          * **self.ott_mrca**: OToL taxon Id for the most recent common ancestor of the ingroup
          * **self.orig_seqlen**: list of the original sequence length of the input data
          * **self.gi_dict**: dictionary, that has all information from sequences found during the blasting.
            * key: GenBank sequence identifier
            * value: dictionary, content depends on blast option, differs between webquery and local blast queries
                * keys - value pairs for local blast:
                    * '^ncbi:gi': GenBank sequence identifier
                    * 'accession': GenBank accession number
                    * 'staxids': Taxon identifier
                    * 'sscinames': Taxon species name
                    * 'pident': Blast  percentage of identical matches
                    * 'evalue': Blast e-value
                    * 'bitscore': Blast bitscore, used for FilterBlast
                    * 'sseq': corresponding sequence
                    * 'title': title of Genbank sequence submission
                * key - values for web-query:
                    * 'accession':Genbank accession number
                    * 'length': length of sequence
                    * 'title': string combination of hit_id and hit_def
                    * 'hit_id': string combination of gi id and accession number
                    * 'hsps': Bio.Blast.Record.HSP object
                    * 'hit_def': title from GenBank sequence
                * optional key - value pairs for unpublished option:
                    * 'localID': local sequence identifier
          * **self._reconciled**: True/False,
          * **self.unpubl_otu_json**: optional, will contain the OTU-dict for unpublished data, if that option is used

        Following functions are called during the init-process:

          * **self._reconcile()**:
                removes taxa, that are not found in both, the phylogeny and the aln
          * **self._reconcile_names()**:
                is used for the own file stuff, it removes the character 'n' from tip names that start with a number

        The physcraper class is then updating:

          * self.aln, self.tre and self.otu_dict, self.ps_otu, self.gi_dict
    """

    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, config_obj,
                 schema=None, taxon_namespace=None):
        debug("build ATT class")
        self.aln = alignment
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix), \
                ("your aln '{}' is not a DnaCharacterMatrix".format(alignment))
        if schema is None:
            self.tre = Tree.get(data=newick,
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        else:
            self.tre = Tree.get(data=newick,
                                schema=schema,
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        assert (self.tre.taxon_namespace is self.aln.taxon_namespace), "tre and aln taxon_namespace are not identical"
        assert isinstance(otu_dict, dict), ("otu_dict '{}' is not of type dict".format(otu_dict))
        self.otu_dict = otu_dict
        self.config = config_obj
        self.ps_otu = 1  # iterator for new otu IDs
        self._reconcile()
        self._reconcile_names()
        self.workdir = os.path.abspath(workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca), ("your ingroup_mrca '{}' is not an integer.".format(ingroup_mrca))
        self.ott_mrca = ingroup_mrca  # ott_ingroup mrca can be pulled directly from phylesystem
        self.orig_seqlen = []  # will get filled in later...
        self.gb_dict = {}  # has all info about new blast seq
        self._reconciled = False
        self.unpubl_otu_json = None

    def _reconcile(self):
        """Taxa that are only found in the tree, or only in the alignment are deleted.

        This checks that the tree "original labels" from phylesystem
        align with those found in the alignment.
        """
        debug("reconcile")
        treed_tax = set()
        for leaf in self.tre.leaf_nodes():
            treed_tax.add(leaf.taxon)
        aln_tax = set()
        for tax, seq in self.aln.items():
            aln_tax.add(tax)
        prune = treed_tax ^ aln_tax
        missing = [i.label for i in prune]
        if missing:
            errmf = 'NAME RECONCILIATION Some of the taxa in the tree are not in the alignment or vice versa' \
                    ' and will be pruned. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            sys.stderr.write(errm)
        del_aln = []
        del_tre = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_tax)
            if taxon in aln_tax:
                # debug(taxon)
                del_aln.append(taxon)
            if taxon in treed_tax:
                del_tre.append(taxon)
        # debug(del_aln)
        # debug(del_tre)
        self.aln.remove_sequences(del_aln)
        self.tre.prune_taxa(del_tre)
        for tax in prune:
            # potentially slow at large number of taxa and large numbers to be pruned
            found = 0
            for otu in self.otu_dict:
                # debug([otu, tax.label])
                if "^ot:originalLabel" in self.otu_dict[otu]:
                    if self.otu_dict[otu][u'^ot:originalLabel'] == tax.label:
                        self.otu_dict[otu]['^physcraper:status'] = "deleted in reconciliation"
                        found = 1
                elif otu == tax.label:
                    self.otu_dict[otu]['^physcraper:status'] = "deleted in reconciliation"
                    found = 1
            if found == 0:
                sys.stderr.write("lost taxon {} in reconcilliation \n".format(tax.label))
            self.aln.taxon_namespace.remove_taxon(tax)
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def _reconcile_names(self):
        """It rewrites some tip names, which kept being an issue when it starts with a number at the beginning.
        Then somehow a n was added to the tip names.

        :return: replaced tip names
        """
        for tax in self.aln.taxon_namespace:
            if tax.label in self.otu_dict.keys():
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for otu in self.otu_dict:
                    original = self.otu_dict[otu].get("^ot:originalLabel")
                    if original == tax.label or original == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))

    def prune_short(self):
        """Prunes sequences from alignment if they are shorter than specified in the config file,
         or if tip is only present in tre.

        Sometimes in the de-concatenating of the original alignment taxa with no sequence are generated
        or in general if certain sequences are really short. This removes those from both the tre and the alignment.

        has test: test_prune_short.py

        :return: prunes aln and tre
        """
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in
                            self.aln]
        avg_seqlen = sum(self.orig_seqlen) / len(self.orig_seqlen)
        #debug("average sequence length is {}".format(avg_seqlen))
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        prune = []
        aln_ids = set()
        for tax, seq in self.aln.items():
            aln_ids.add(tax.label)
            if len(seq.symbols_as_string().translate(None, "-?")) <= seq_len_cutoff:
                prune.append(tax)
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short "
                     "step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                self.remove_taxa_aln_tre(tax.label)
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]["^physcraper:status"] = "deleted in prune short"

        # self.trim()
        # out-commented next assert line, as this does not run if we prune aln before placing new seq in tre
        # debug(self.aln.taxon_namespace)
        # debug(self.tre.taxon_namespace)
        # assert self.aln.taxon_namespace == self.tre.taxon_namespace
        # next function has assert statement included
        self.check_tre_in_aln()
        self._reconciled = 1

    def trim(self):
        """ It removes bases at the start and end of alignments, if they are represented by less than the value
        specified in the config file. E.g. 0.75 given in config means, that 75% of the sequences need to have a
        base present

        Ensures, that not whole chromosomes get dragged in. It's cutting the ends of long sequences.

        has test: test_trim.py
        """
        # debug('in trim')
        taxon_missingness = self.config.trim_perc
        i = 0
        seqlen = len(self.aln[i])
        while seqlen == 0:
            i = i + 1
            seqlen = len(self.aln[i])
        for tax in self.aln:
            if len(self.aln[tax]) != seqlen:
                sys.stderr.write("can't trim un-aligned inputs, moving on")
                return
        start = 0
        stop = seqlen
        cutoff = len(self.aln) * taxon_missingness
        for i in range(seqlen):
            counts = {"?": 0, "-": 0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ["?", "-"]:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:  # first ok column
                start = i
                break
        for i in range(seqlen, 0, -1):  # previously seqlen-1, that cuts off last character of aln, I changed it.
            counts = {'?': 0, '-': 0}
            for tax in self.aln:
                call = self.aln[tax][i - 1].label  # changing seqlen-1 to seqlen requires that we have here i-1
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                stop = i
                break
        # here alignment gets shortened to start:stop
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        # make sure that tre is presented in aln
        self.check_tre_in_aln()
        if _VERBOSE:
            sys.stdout.write("trimmed alignment ends to < {} missing taxa, "
                             "start {}, stop {}\n".format(taxon_missingness, start, stop))
        return

    def check_tre_in_aln(self):
        """Makes sure that everything which is in tre is also found in aln.

        Extracted method from trim. Not sure we actually need it there.
        """
        aln_ids = set()
        for taxon in self.aln:
            aln_ids.add(taxon.label)
        assert aln_ids.issubset(self.otu_dict.keys()), \
            ([x for x in aln_ids if x not in self.otu_dict.keys()], self.otu_dict.keys())

        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        for leaf in self.tre.leaf_nodes():
            if leaf.taxon not in aln_ids:
                debug(leaf.taxon)
                self.tre.prune_taxa([leaf])
                self.tre.prune_taxa_with_labels([leaf.taxon])
                self.tre.prune_taxa_with_labels([leaf])
                treed_taxa.remove(leaf.taxon)
        assert treed_taxa.issubset(aln_ids), (treed_taxa, aln_ids)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln and tre and updates otu_dict,
        takes a single taxon_label as input.

        note: has test, test_remove_taxa_aln_tre.py

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return: removes information/data from taxon_label
        """
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.discard_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found
            # the first prune does not remove it sometimes...
        if tax2:
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, updated otu_dict but was never in tre or aln!"

    def add_otu(self, gb_id, ids_obj):
        """ Generates an otu_id for new sequences and adds them into self.otu_dict.
        Needs to be passed an IdDict to do the mapping.

        :param gb_id: the Genbank identifier/ or local unpublished
        :param ids_obj: needs to IDs class to have access to the taxonomic information
        :return: the unique otu_id - the key from self.otu_dict of the corresponding sequence
        """
        # debug("add_otu function")
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        ott_id = None
     #   debug("trying to add an otu with accesion {}".format(gb_id))
        ncbi_id, tax_name = get_tax_info_from_acc(gb_id, self, ids_obj)
     #   debug("ADD OTU: accession {} ncbi_id {}, taxon_name {}".format(gb_id, ncbi_id, tax_name))
        if ncbi_id == None:
            return None
        if ncbi_id in ids_obj.ncbi_to_ott.keys():
            ott_id = int(ids_obj.ncbi_to_ott[int(ncbi_id)])
        else:
            #debug("{} Ncbi id not found in ott_ncbi dictionaries\n".format(ncbi_id))
            ott_id = None
        if ott_id in ids_obj.ott_to_name:
            ott_name = ids_obj.ott_to_name[ott_id]
        else:
            ott_name = None
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]["^ncbi:title"] = self.gb_dict[gb_id]["title"]
        self.otu_dict[otu_id]["^ncbi:taxon"] = ncbi_id
        self.otu_dict[otu_id]["^ncbi:TaxonName"] = tax_name
        self.otu_dict[otu_id]["^ot:ottId"] = ott_id
        self.otu_dict[otu_id]["^physcraper:status"] = "query"
        self.otu_dict[otu_id]["^ot:ottTaxonName"] = ott_name
        # last_blasted date infos: 1800 = never blasted; 1900 = blasted 1x, not added; this century = blasted and added
        self.otu_dict[otu_id]["^physcraper:last_blasted"] = "1800/01/01"
        if gb_id[:6] == "unpubl":
            self.otu_dict[otu_id]["^physcraper:status"] = "local seq"
            self.otu_dict[otu_id]["^ot:originalLabel"] = self.gb_dict[gb_id]["localID"]
            self.otu_dict[otu_id]['^user:TaxonName'] = self.gb_dict[gb_id][u'^user:TaxonName']
        else:
            self.otu_dict[otu_id]["^ncbi:gi"] = self.gb_dict[gb_id]["^ncbi:gi"]
            self.otu_dict[otu_id]["^ncbi:accession"] = gb_id
        # get a name for the OTU, no matter from which source
        if tax_name is not None:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = tax_name
        elif ott_name is not None:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = ott_name
        elif self.otu_dict[otu_id].get('^user:TaxonName'):
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = self.otu_dict[otu_id]['^user:TaxonName']
        else:
            self.otu_dict[otu_id]["^physcraper:TaxonName"] = "ACC_{}".format(gb_id)
        assert self.otu_dict[otu_id]["^physcraper:TaxonName"]  # is not None
        if _DEBUG >= 2:
            sys.stderr.write("acc:{} assigned new otu: {}\n".format(gb_id, otu_id))
        return otu_id

    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """This writes out needed files for papara (except query sequences).
        Papara is finicky about trees and needs phylip format for the alignment.

        NOTE: names for tree and aln files should not be changed, as they are hardcoded in align_query_seqs().

        Is only used within func align_query_seqs.
        """
        #debug('write papara files')
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")  # Papara is diffffffficult about root
        tmptre = tmptre.replace("'", "_")
        fi = open("{}/{}".format(self.workdir, treefilename), "w")
        fi.write(tmptre)
        fi.close()
        self.aln.write(path="{}/{}".format(self.workdir, alnfilename), schema="phylip")

    def write_files(self, treepath="physcraper.tre", treeschema="newick", alnpath="physcraper.fas", alnschema="fasta"):
        """Outputs both the streaming files, labeled with OTU ids.
        Can be mapped to original labels using otu_dict.json or otu_seq_info.csv"""
        #debug("write_files")
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema, unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)

    def write_labelled(self, label, treepath=None, alnpath=None, norepeats=True, add_gb_id=False):
        """output tree and alignment with human readable labels
        Jumps through a bunch of hoops to make labels unique.

        NOT MEMORY EFFICIENT AT ALL

        Has different options available for different desired outputs

        :param label: which information shall be displayed in labelled files: possible options:
                    '^ot:ottTaxonName', '^user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"
        :param treepath: optional: full file name (including path) for phylogeny
        :param alnpath:  optional: full file name (including path) for alignment
        :param norepeats: optional: if there shall be no duplicate names in the labelled output files
        :param add_gb_id: optional, to supplement tiplabel with corresponding GenBank sequence identifier
        :return: writes out labelled phylogeny and alignment to file
        """
        #debug("write labelled files")
        if treepath is None:
            treepath = "{}/{}".format(self.workdir, "labelled.tre")
        if alnpath is None:
            alnpath = "{}/{}".format(self.workdir, 'labelled.aln')
        assert label in ['^ot:ottTaxonName', '^user:TaxonName', '^physcraper:TaxonName',
                         "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        tmp_fasta = self.aln.as_string(schema="fasta")
        tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta,
                                         schema="fasta",
                                         taxon_namespace=tmp_tre.taxon_namespace)
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label is None:
                if self.otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_{}".format(self.otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"),
                                                            self.otu_dict[taxon.label].get('^physcraper:TaxonName', "unk"))
            new_label = str(new_label).replace(' ', '_')
            if add_gb_id:
                gb_id = self.otu_dict[taxon.label].get('^ncbi:accession')
                if gb_id is None:
                    gb_id = self.otu_dict[taxon.label].get("^ot:originalLabel")
                new_label = "_".join([new_label, str(gb_id)])
                sp_counter = 2
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, str(sp_counter)])
                    sp_counter += 1
            else:
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, taxon.label])
            taxon.label = new_label
            new_names.add(new_label)
        tmp_tre.write(path=treepath,
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path=alnpath,
                      schema="fasta")

    def write_otus(self, filename, schema="table"):
        """Writes out OTU dict as json or table.

        :param filename: filename
        :param schema: either table or json format
        :return: writes out otu_dict to file
        """
        # TODO: schema is unused!
        assert schema in ["table", "json"]
        with open("{}/{}".format(self.workdir, filename), "w") as outfile:
            json.dump(self.otu_dict, outfile)

    def dump(self, filename=None):
        """writes pickled files from att class"""
        if filename:
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/att_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi)


#####################################
def get_nexson(study_id, phylesystem_loc):
    """Grabs nexson from phylesystem"""
    phy = PhylesystemAPI(get_from=phylesystem_loc)
    nexson = phy.get_study(study_id)["data"]
    return nexson


def get_mrca_ott(ott_ids):
    """finds the mrca of the taxa in the ingroup of the original
    tree. The blast search later is limited to descendants of this
    mrca according to the ncbi taxonomy

    Only used in the functions that generate the ATT object.

    :param ott_ids: list of all OToL identifiers for tip labels in phylogeny
    :return: OToL identifier of most recent common ancestor or ott_ids
    """
    debug("get_mrca_ott")
    if None in ott_ids:
        ott_ids.remove(None)
    synth_tree_ott_ids = []
    ott_ids_not_in_synth = []
    for ott in ott_ids:
        r = opentree_helpers.check_if_ottid_in_synth(ott)
        if r == 1:
            synth_tree_ott_ids.append(ott)
        else:
            ott_ids_not_in_synth.append(ott)
    if len(synth_tree_ott_ids) == 0:
        sys.stderr.write('No sampled taxa were found in the current synthetic tree. '
                         'Please find and input and appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit(-3)
    mrca_node = tree_of_life.mrca(ott_ids=synth_tree_ott_ids, wrap_response=False)  # need to fix wrap eventually
    if u'nearest_taxon' in mrca_node.keys():
        tax_id = mrca_node[u'nearest_taxon'].get(u'ott_id')
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node[u'nearest_taxon'][u'name']))
    elif u'taxon' in mrca_node['mrca'].keys():
        tax_id = mrca_node['mrca'][u'taxon'][u'ott_id']
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node['mrca'][u'taxon'][u'name']))
    else:
        # debug(mrca_node.keys())
        sys.stderr.write('(v3) MRCA of sampled taxa not found. Please find and input an '
                         'appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit(-4)
    return tax_id

#####################################


class IdDicts(object):
    """Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

        To build the class the following is needed:

          * **config_obj**: Object of class config (see above)
          * **workdir**: the path to the assigned working directory

        During the initializing process the following self objects are generated:

          * **self.workdir**: contains path of working directory
          * **self.config**: contains the Config class object
          * **self.ott_to_ncbi**: dictionary

              * key: OToL taxon identifier
              * value: ncbi taxon identifier
          * **self.ncbi_to_ott**: dictionary

              * key: OToL taxon identifier
              * value: ncbi taxon identifier
          * **self.ott_to_name**: dictionary

              * key: OToL taxon identifier
              * value: OToL taxon name
          * **self.acc_ncbi_dict**: dictionary

              * key: Genbank identifier
              * value: ncbi taxon identifier
          * **self.spn_to_ncbiid**: dictionary

              * key: OToL taxon name
              * value: ncbi taxon identifier
          * **self.ncbiid_to_spn**: dictionary

              * key: ncbi taxon identifier
              * value: ncbi taxon name

          user defined list of mrca OTT-ID's #TODO this is flipped form the dat aobj .ott_mrca. On purpose?
         #reomved mrca's from ida, and put them into scrape object

          * **Optional**:

              * depending on blasting method:
               * self.ncbi_parser: for local blast, initializes the ncbi_parser class, that contains information about rank and identifiers
               * self.otu_rank: for remote blast to store the rank information
    """

    def __init__(self, config_obj, workdir):
        """Generates a series of name disambiguation dicts"""
        self.config = config_obj
        assert self.config.email
        self.ott_to_ncbi = {} 
        self.ncbi_to_ott = {}  # used to get ott_id for new Genbank query taxa
        self.ott_to_name = {}  # used in add_otu to get name from otuId
        self.acc_ncbi_dict = {}  # filled by ncbi_parser (by subprocess in earlier versions of the code).
        self.spn_to_ncbiid = {}  # spn to ncbi_id, it's only fed by the ncbi_data_parser, but makes it faster
        self.ncbiid_to_spn = {} #TODO when is this generated? MK: well, here. it is filled with information from genbank to speed up translation between ncbi_taxon_ids and names. similar to  acc_ncbi_dict and spn_to_ncbiid.

        fi = open(config_obj.ott_ncbi)  # This is in the taxonomy folder of the repo, needs to be updated by devs when OpenTree taxonomy changes.
        for lin in fi:  # TODO This is insanely memory inefficient, how about using a pandas dataframe?
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1]) #TODO EJM if an ottid matches 2 ncbi ids, only 1 appears
            self.ncbi_to_ott[int(lii[1])] = int(lii[0]) #TODO EJM - ott ids that match two ncbi ids just show up once
            self.ott_to_name[int(lii[0])] = lii[2].strip()  # todo merge into ott_to_ncbi?
            assert len(self.ott_to_ncbi) > 0
            assert len(self.ncbi_to_ott) > 0
            assert len(self.ott_to_name) > 0
        fi.close()
        # TODO: pandas solution? requires to rewrite usages of self.ott_to_ncbi, self.ncbi_to_ott, self.ott_to_name
        if config_obj.blast_loc == 'remote':
            self.otu_rank = {}  # used only for web queries - contains taxonomic hierarchy information
        else:  # ncbi parser contains information about spn, tax_id, and ranks
            self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                                       nodes_file=self.config.ncbi_parser_nodes_fn)
        self.acc_tax_seq_dict = {} 

#moved self.get_ncbi_mrca to physcraper data object


# removed funtion ottid_to_ncbiid.
# if ott is isn't mapped to ncbi_id in the opentree taxonomy, (ott_to_ncbi) we shouldn't map it here.

    def get_ncbiid_from_tax_name(self, tax_name): 
    #TODO when do we need this? Is this just when we read in the original data?
        """Get the ncbi_id from the species name using ncbi web query.

        :param tax_name: species name
        :return: corresponding ncbi id
        """
        ncbi_id = None
        if tax_name in self.spn_to_ncbiid:
            ncbi_id = self.spn_to_ncbiid[tax_name]
        else:
            try:
                tries = 15
                for i in range(tries):
                    try:
                        Entrez.email = self.config.email
                        if tries >= 5:
                            ncbi_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
                        else:
                            tax_name = "'{}'".format(tax_name)
                            ncbi_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]

                        ncbi_id = int(ncbi_id)
                    except (IndexError, HTTPError) as err:
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            raise
                    break
            except (IndexError, HTTPError) as err:
                #debug("except")
                try:
                    ncbi = NCBITaxa()
                    tax_info = ncbi.get_name_translator([tax_name])
                    #debug(tax_info)
                    if tax_info == {}:
                        tax_name = "'{}'".format(tax_name)
                        tax_info = ncbi.get_name_translator([tax_name])
                    ncbi_id = int(tax_info.items()[0][1][0])
                except (IndexError, HTTPError) as err:
                    sys.stderr.write("Taxon name does not match any name in ncbi. Check that name is written "
                                     "correctly: {}! We set it to unidentified".format(tax_name))
                    tax_name = 'unidentified'
                    ncbi_id = 0
        assert type(ncbi_id) is int
        self.spn_to_ncbiid[tax_name] = ncbi_id
        return ncbi_id

    def get_rank_info_from_web(self, taxon_name):
        #TODO, why input name rather than ID here?
        """Collects rank and lineage information from ncbi,
        used to delimit the sequences from blast,
        when the web blast service is used.
        """
        tax_name = taxon_name.replace(" ", "_")
        # if tax_name not in self.otu_rank.keys():
        ncbi_id = self.get_ncbiid_from_tax_name(tax_name)
        if ncbi_id == None:
            self.otu_rank[ncbi_id] = {"taxon id": ncbi_id, "lineage": 'life', "rank": 'unassigned'}
        else:
            ncbi = NCBITaxa()
            lineage = ncbi.get_lineage(ncbi_id)
            lineage2ranks = ncbi.get_rank(lineage)
            tax_name = str(tax_name).replace(" ", "_")
            assert type(ncbi_id) is int
            self.otu_rank[ncbi_id] = \
                {"taxon id": ncbi_id, "lineage": lineage, "rank": lineage2ranks, "taxon name": tax_name}
        return ncbi_id

 #removed function find_tax_id because it wasn't being used
   
    def get_tax_seq_acc(self, acc=None):
        gb_id = acc
        if len(gb_id.split(".")) == 1:
            debug("accession number {} not recognized".format(gb_id))
        if gb_id in self.acc_tax_seq_dict:
            tax_name = self.acc_tax_seq_dict[gb_id]["taxname"]
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:taxon"]
            seq = self.acc_tax_seq_dict[gb_id]["seq"]
        else:
            read_handle = self.entrez_efetch(gb_id)
            tax_name = get_ncbi_tax_name(read_handle)
            ncbi_id = get_ncbi_tax_id(read_handle)
            seq = read_handle[0][u'GBSeq_sequence']
            tax_name = tax_name.replace(" ","_") #TODO check that searches are using names without spaces 
            self.ncbiid_to_spn[ncbi_id] = tax_name
            self.acc_ncbi_dict[gb_id] = ncbi_id
            self.acc_tax_seq_dict[gb_id] = {'taxname':tax_name, "^ncbi:taxon":ncbi_id, 'seq':seq} #This is going to be a memory hog...
        assert ncbi_id is not None
        return ncbi_id, tax_name, seq

    def get_ncbiid_from_acc(self, acc):
        '''checks local dicts, and then runs eftech to get ncbi id for accession'''
        gb_id = acc
        if gb_id in self.acc_ncbi_dict:#TODO if the accession number and tax id are here, does that mean the name is in ncbiid_to_spn?
            ncbi_id = self.acc_ncbi_dict[gb_id]
        elif gb_id in self.acc_tax_seq_dict:
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:taxon"]
        else:
            taxid, taxname, seq = self.get_tax_seq_acc(self, acc=gb_id)
            ncbi_id = taxid
        return ncbi_id

   
# If we didn't find name when creating otu_dict entry, we shouldn't look again.
# if we have accession number, we can find name using get_tax_info_from_acc
    def find_name_acc(self, acc):
        ncbi_id = self.get_ncbiid_from_acc(acc)
        tax_name = self.ncbiid_to_spn[ncbi_id]  #TODO check that values are getting enetered into this dict
        return tax_name

    def find_name_otu(self, otu_dict_entry=None):
        """ Find the taxon name in the  otu_dict entry or of a Genbank accession number.
        If not already known it will ask ncbi using the accession number

        :param otu_dict_entry: otu_dict entry
        :param acc: Genbank accession number
        :return: ncbi taxon name
        """
        # debug("find_name")
        inputinfo = False
        if otu_dict_entry is not None:
            inputinfo = True
        assert inputinfo is True
        tax_name = None
        ncbi_id = None
        if otu_dict_entry:
            # debug(otu_dict_entry)
            if "^physcraper:TaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^physcraper:TaxonName"]
            elif "^ot:ottTaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^ot:ottTaxonName"]
            elif "^user:TaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^user:TaxonName"]
        assert tax_name is not None
        tax_name = tax_name.replace(" ", "_")
        return tax_name


    def entrez_efetch(self, gb_id):
        """ Wrapper function around efetch from ncbi to get taxonomic information if everything else is failing.
            Also used when the local blast files have redundant information to access the taxon info of those sequences.
        It adds information to various id_dicts.

        :param gb_id: Genbank identifier
        :return: read_handle
        """
        tries = 10
        Entrez.email = self.config.email
        handle = None

        # method needs delay because of ncbi settings
        for i in range(tries):
            try:
                # print("try")
                delay = 1.0
                previous = time.time()
                while True:
                    current = time.time()
                    wait = previous + delay - current
                    if wait > 0:
                        # print("if", wait)
                        time.sleep(wait)
                        previous = current + wait
                    else:
                        # print("else", wait)
                        previous = current
                    if delay + .5 * delay <= 5:
                        # print("if2", delay)
                        delay += .5 * delay
                    else:
                        # print("else2",  delay)
                        delay = 5
                    # print("read handle")
                    handle = Entrez.efetch(db="nucleotide", id=gb_id, retmode="xml")
                    assert handle is not None, ("your handle file to access data from efetch does not exist. "
                                                "Likely an issue with the internet connection of ncbi. Try rerun...")
                    read_handle = Entrez.read(handle)
                    handle.close()

                    return read_handle
            except (IndexError, HTTPError) as e:
                if i < tries - 1:  # i is zero indexed
                    continue
                else:
                    raise
            # break
        assert handle is not None, ("your handle file to access data from efetch does not exist. "
                                    "Likely an issue with the internet connection of ncbi. Try rerun...")
        read_handle = Entrez.read(handle)
        handle.close()
        return read_handle


    def dump(self, filename=None):
        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("id_pickle.p", "wb")
        pickle.dump(self, ofi)


class PhyscraperScrape(object):
    """
    This is the class that does the perpetual updating

        To build the class the following is needed:

          * **data_obj**: Object of class ATT (see above)
          * **ids_obj**: Object of class IdDict (see above)

        During the initializing process the following self.objects are generated:

          * **self.workdir**: path to working directory retrieved from ATT object = data_obj.workdir
          * **self.logfile**: path of logfile
          * **self.data**: ATT object
          * **self.ids**: IdDict object
          * **self.config**: Config object
          * **self.new_seqs**: dictionary that contains the newly found seq using blast:

            * key: gi id
            * value: corresponding seq
          * **self.new_seqs_otu_id**: dictionary that contains the new sequences that passed the remove_identical_seq() step:

            * key: otu_id
            * value: see otu_dict, is a subset of the otu_dict, all sequences that will be newly added to aln and tre
          * **self.otu_by_gi**: dictionary that contains ????:

            * key:
            * value:
          * **self._to_be_pruned**: list that contains ????
          * **self.mrca_ott_list**: user defined list of mrca OTT-ID's
          * **self.mrca_ncbi_list**:  int or list of ncbi identifier of mrca

          * **self.tmpfi**: path to a file or folder???
          * **self.blast_subdir**: path to folder that contains the files writen during blast

          * **self.newseqs_file**: filename of files that contains the sequences from self.new_seqs_otu_id
          * **self.date**: Date of the run - may lag behind real date!
          * **self.repeat**: either 1 or 0, it is used to determine if we continue updating the tree, no new seqs found = 0
          * **self.newseqs_acc**: list of all gi_ids that were passed into remove_identical_seq(). Used to speed up adding process
          * **self.blacklist**: list of gi_id of sequences that shall not be added or need to be removed. Supplied by user.
          * **self.seq_filter**: list of words that may occur in otu_dict.status and which shall not be used in the building of FilterBlast.sp_d (that's the main function), but it is also used as assert statement to make sure unwanted seqs are not added.
          * **self.unpublished**: True/False. Used to look for local unpublished seq that shall be added if True.
          * **self.path_to_local_seq:** Usually False, contains path to unpublished sequences if option is used.

        Following functions are called during the init-process:

            * **self.reset_markers()**: adds things to self: I think they are used to make sure certain function run, if program crashed and pickle file is read in.
                * self._blasted: 0/1, if run_blast_wrapper() was called, it is set to 1 for the round.
                * self._blast_read: 0/1, if read_blast_wrapper() was called, it is set to 1 for the round.
                * self._identical_removed: 0
                * self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
                * self._query_seqs_aligned: 0
                * self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
                * self._reconciled: 0
                * self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
            * **self.OToL_unmapped_tips()**: function that either removes or maps unmapped taxa from OToL studies
    """
    def __init__(self, data_obj, ids_obj, ingroup_mrca=None):
        assert isinstance(data_obj, AlignTreeTax)
        assert isinstance(ids_obj, IdDicts)
        debug("start base class init")
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
#        assert data_obj.config == ids_obj.config
        self.config = self.ids.config  # pointer to config
        self.new_seqs = {}  # all new seq after read_blast_wrapper
        self.new_seqs_otu_id = {}  # only new seq which passed remove_identical
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)  # TODO: For what do we want to use this? Unused!
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = ""
        self.date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        self.repeat = 1  # used to determine if we continue updating the tree
        self.newseqs_acc = []  # all ever added Genbank accession numbers during any PhyScraper run, used to speed up adding process
        self.blacklist = []  # remove sequences by default
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,",
                           "local"]  # TODO MK: try to move completely to FilterBlast class
        self.reset_markers()
        self.unpublished = False  # used to look for local unpublished seq that shall be added.
        self.path_to_local_seq = False  # path to unpublished seq.
        self.backbone = False
        self.OToL_unmapped_tips()  # added to do stuff with un-mapped tips from OToL
        self.ids.ingroup_mrca = data_obj.ott_mrca  # added for mrca ingroup list
        self.gb_not_added = []  # list of blast seqs not added
        self.del_superseq = set()  # items that were deleted bc they are superseqs, needed for assert statement
        self.mrca_ott_list = ingroup_mrca
        assert type(self.mrca_ott_list) in [int, list] or self.mrca_ott_list is None
        self.mrca_ncbi_list = set()  # corresponding ids for mrca_ott list
        if self.mrca_ott_list is not None:
            self.get_ncbi_mrca()
        # self.mrca_ncbi_of_all = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        if len(self.mrca_ncbi_list) == 0:
            # self.mrca_ncbi = set()
            self.mrca_ncbi_list.add(ids_obj.ott_to_ncbi[data_obj.ott_mrca])
        #debug(self.mrca_ncbi_list)
        debug("created physcaper ncbi_mrca {}, len {}".format(self.mrca_ncbi_list, len(self.mrca_ncbi_list)))
        self.map_taxa_to_ncbi()


    def get_ncbi_mrca(self):
        """ get the ncbi tax ids from a set of mrca ott ids.
        """
        if type(self.mrca_ott_list) is not int:
            for ott_id in self.mrca_ott_list:
                ncbi_id = self.ottid_to_ncbiid(ott_id)
                if ncbi_id is not None:
                    self.mrca_ncbi_list.add(ncbi_id)
                else:
                    sys.stderr.write("ott_id {} is not mapped to ncbi_id".format(ott_id))
        else:
            ncbi_id = self.ottid_to_ncbiid(self.mrca_ott_list)
            if ncbi_id is not None:
                self.mrca_ncbi_list.add(ncbi_id)
                
    def map_taxa_to_ncbi(self):
        for otu in self.data.otu_dict:
            if self.data.otu_dict[otu].get("^ncbi:taxon") == None:
                if self.data.otu_dict[otu].get("^ot:ottId"):
                    ottid = self.data.otu_dict[otu]["^ot:ottId"]
                    self.data.otu_dict[otu]["^ncbi:taxon"]=self.ids.ott_to_ncbi.get(ottid,0)


    def ottid_to_ncbiid(self, ott_id):
        """ Find ncbi id for ott id.
        Is only used for the mrca list thing!
        """
        debug("ottid to ncbiid")
        #debug(ott_id)
        ncbi_id = None
        if ott_id in self.ids.ott_to_ncbi:
            ncbi_id = self.ids.ott_to_ncbi[ott_id]
        else:
            sys.stderr.write("ottid {} does not have a direct match to ncbi.".format(ott_id))
            assert ott_id not in self.ids.ott_to_name
            #if ott_id is not in ott_to_ncbi, it doean't have a direct map to ncbi
            #and it shouldn't be in ott_name
            #
        return ncbi_id



    # TODO is this the right place for this? MK: According to PEP8, no...
    def reset_markers(self):
        self._blasted = 0
        self._blast_read = 0
        self._query_seqs_written = 0
        self._query_seqs_placed = 0
        self._full_tree_est = 0

    def reset_new_seqs_acc(self):
        """ Needs to be reseted if you want to rerun the filtering to get lower rank taxa added"""
        self.newseqs_acc = []

    def OToL_unmapped_tips(self):
        """Assign names or remove tips from aln and tre that were not mapped during initiation of ATT class.
        """
        debug("OTOL unmapped")
        if self.config.unmapped == "remove":
            for key in self.data.otu_dict:
                if "^ot:ottId" not in self.data.otu_dict[key]:
                    # second condition for OToL unmapped taxa, not present in own_data
                    if u"^ot:treebaseOTUId" in self.data.otu_dict[key]:
                        self.data.remove_taxa_aln_tre(key)
        elif self.config.unmapped == "keep":
            i = 1
            for key in self.data.otu_dict:
                i = i + 1
                if "^ot:ottId" not in self.data.otu_dict[key]:
                    self.data.otu_dict[key]["^ot:ottId"] = self.data.ott_mrca
                    if self.data.ott_mrca in self.ids.ott_to_name:
                        self.data.otu_dict[key]['^ot:ottTaxonName'] = self.ids.ott_to_name[self.data.ott_mrca]
                    else:
                        debug("think about a way...")
                        tx = APIWrapper().taxomachine
                        nms = tx.taxon(self.data.ott_mrca)
                        taxon_name = nms[u'unique_name']
                        self.data.otu_dict[key]['^ot:ottTaxonName'] = "unknown_{}".format(taxon_name)

    def run_local_blast_cmd(self, query, taxon_label, fn_path):
        """Contains the cmds used to run a local blast query, which is different from the web-queries.

        :param query: query sequence
        :param taxon_label: corresponding taxon name for query sequence
        :param fn_path: path to output file for blast query result

        :return: runs local blast query and writes it to file
        """
        abs_blastdir = os.path.abspath(self.blast_subdir)
        abs_fn = os.path.abspath(fn_path)
        toblast = open("{}/tmp.fas".format(os.path.abspath(self.blast_subdir)), "w+")
        toblast.write(">{}\n".format(taxon_label))
        toblast.write("{}\n".format(query))
        toblast.close()
        assert os.path.isdir(self.config.blastdb), ("blast dir does not exist: '{}'.".format(self.config.blastdb))
        with cd(self.config.blastdb):
            # this format (6) allows to get the taxonomic information at the same time
            outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid'"
            # outfmt = " -outfmt 5"  # format for xml file type
            # TODO query via stdin
            blastcmd = "blastn -query " + "{}/tmp.fas".format(abs_blastdir) + \
                       " -db {}/nt -out ".format(self.config.blastdb) + abs_fn + \
                       " {} -num_threads {}".format(outfmt, self.config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(self.config.hitlist_size,
                                                                  self.config.hitlist_size)
            os.system(blastcmd)

    def local_blast_for_unpublished(self, query, taxon):
        """
        Run a local blast search if the data is unpublished.

        :param query: query sequence
        :param taxon: taxon.label used as identifier for the sequences
        :return: xml files with the results of the local blast
        """
        with cd(os.path.join(self.workdir, "blast")):
            debug("run against local unpublished data")
            debug(self.blast_subdir)
            toblast = open("{}/tmp.fas".format(self.blast_subdir), "w")
            toblast.write(">{}\n".format(taxon))
            toblast.write("{}\n".format(query))
            toblast.close()
            blast_db = "local_unpubl_seq_db"
            output = "tst_fn"
            blastcmd = "blastn -query {}/tmp.fas -db {} -out output_{}.xml " \
                       "-outfmt 5".format(self.blast_subdir, blast_db, output)
            os.system(blastcmd)

    def run_web_blast_query(self, query, equery, fn_path):
        """Equivalent to run_local_blast_cmd() but for webqueries,
        that need to be implemented differently.

        :param query: query sequence
        :param equery: method to limit blast query to mrca
        :param fn_path: path to output file for blast query result
        :return: runs web blast query and writes it to file
        """
        if self.config.url_base:
            debug("blasting {} using {}".format(self.config.url_base))
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          url_base=self.config.url_base,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size,
                                          num_threads=self.config.num_threads)
        else:
            debug("blasting {} using webservice".format(fn_path))
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size)
        save_file = open(fn_path, "w")
        save_file.write(result_handle.read())
        result_handle.close()
        save_file.close()

    def run_blast_wrapper(self):  # TODO Should this be happening elsewhere?
        """generates the blast queries and saves them depending on the blasting method to different file formats

        It runs blast if the sequences was not blasted since the user defined threshold in the config file (delay).

        :return: writes blast queries to file
        """
        delay = self.config.delay
        debug("run_blast_wrapper")
        debug(self.blast_subdir)
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            if otu_id in self.data.otu_dict:
                if _VERBOSE:
                    sys.stdout.write("blasting {}\n".format(otu_id))
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(
                    last_blast, "%Y/%m/%d")).days)
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                if self.unpublished:
                    self.local_blast_for_unpublished(query, taxon.label)
                    if self.backbone is True:
                        self.data.otu_dict[otu_id]["^physcraper:last_blasted"] = today
                else:
                    if time_passed > delay:
                        if self.config.blast_loc == "local":
                            file_ending = "txt"
                        else:
                            file_ending = "xml"
                        if self.config.gb_id_filename is True:
                            fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label)
                            fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                        else:
                            fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                        # if _DEBUG:
                        #     sys.stdout.write("attempting to write {}\n".format(fn_path))
                        if not os.path.isfile(fn_path):
                            if _VERBOSE:
                                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                            if self.config.blast_loc == 'local':
                                self.run_local_blast_cmd(query, taxon.label, fn_path)
                            if self.config.blast_loc == 'remote':
                                len_ncbi = len(self.mrca_ncbi_list)
                                equery = ''
                                for ncbi_id in self.mrca_ncbi_list:  # add taxids of list to blast search
                                    if len_ncbi >= 2:
                                        equery = equery + "txid{}[orgn] OR ".format(ncbi_id)
                                        len_ncbi = len_ncbi - 1
                                    else:
                                        equery = equery + "txid{}[orgn]) ".format(ncbi_id)
                                equery = "(" + equery + "AND {}:{}[mdat]".format(last_blast, today)
                                self.run_web_blast_query(query, equery, fn_path)
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        else:
                            if _DEBUG:
                                sys.stdout.write("file {} exists in current blast run. Will not blast, "
                                                 "delete file to force\n".format(fn_path))
                            if _DEBUG_MK == 1:
                                self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                    else:
                        if _VERBOSE:
                            sys.stdout.write("otu {} was last blasted {} days ago and is not being re-blasted. "
                                             "Use run_blast_wrapper(delay = 0) to force a search.\n".format(otu_id,
                                                                                                            last_blast))
        self._blasted = 1

    def get_new_seqs_for_mergedseq(self, var_list, query_dict):
        """
        Get the sequence information and fill dict for new sequences when there are merged information.
        
        :param var_list: list of needed variables
        :param query_dict: dict with a query seqs that is filled and returned
        """
        gb_acc, gi_id, sseq, staxids, sscinames, pident, evalue, bitscore, salltitles, sallseqid =  [var for var in var_list]

        found_taxids = set()
        found_spn = set()
        # FOR MERGED SEQS
        staxids_l = staxids.split(";")
        sscinames_l = sscinames.split(";")
        sallseqid_l = sallseqid.split(";")
        salltitles_l = salltitles.split("<>")

        # make sure you have the correct infos
        tax_id_l = self.get_taxid_from_acc(gb_acc)
        # debug([staxids, gb_acc, tax_id_l])
        for item in tax_id_l:
            try:
                assert str(item) in staxids_l, (item, staxids_l, tax_id_l, gb_acc)
            except AssertionError:
                sys.stderr.write("Genbank accession {} does map to different taxonid in database: {} vs {}. Check with ncbi\n".format(gb_acc, item, staxids_l)) 

        # this while loop is here to speed up the process of finding the correct information
        count = 0
        spn_range = 0
        stop_while = False
        id_before = 0
        id_now = -2
        while len(found_taxids) < len(set(staxids_l)):  # as long as i have not found all taxids for the seq
            count += 1
            if stop_while:
                # debug("stop while")
                break
            if count == 5:
                break  # too many tries to find correct number of redundant taxa
            elif count == 1:
                for i in range(0, len(sallseqid_l)):
                    if len(found_taxids) == len(staxids_l):  # if we found all taxon_ids present in the initial list, stop looking for more
                        break
                    gb_acc = get_acc_from_blast(sallseqid_l[i])
                    gi_id = get_gi_from_blast(sallseqid_l[i])
                    stitle = salltitles_l[i]
                    # if both var are the same, we do not need to search GB for taxon info
                    # staxids = tax_id_l[i]
                    staxids = int(self.get_taxid_from_acc(gb_acc)[0])  # corresponding taxid is always in first position
                    # debug(qtaxid)
                    id_before = id_now
                    id_now = staxids
          
                    for item in tax_id_l:
                        try:
                            assert str(staxids) in staxids_l or str(staxids) in tax_id_l, (item, staxids_l, tax_id_l, gb_acc)
                        except AssertionError:
                            sys.stderr.write("Merged Genbank accession {} does map to different taxonid in database: {} vs {}. Check with ncbi\n".format(gb_acc, item, staxids_l)) 
                    # assert str(staxids) in staxids_l or staxids in tax_id_l, (str(staxids), staxids_l, tax_id_l)
                    
                    # if gb acc was already read in before stop the for loop
                    if gb_acc in query_dict or gb_acc in self.data.gb_dict:
                        stop_while = True
                        break
                    sscinames = self.ids.ncbi_parser.get_name_from_id(staxids)
                    # # likely not a speed up to just the single line above
                    # # sometimes if multiple seqs are merged,
                    # # we lack information about which taxon is which gb_acc...
                    # # test it here:
                    # # if we have same number of gb_acc and taxon id go ahead as usual
                    # # debug("get name from id")
                    # if len(sallseqid_l) == len(staxids_l):
                    #     sscinames = sscinames_l[i]
                    # # only one taxon id present, all are from same taxon
                    # elif len(staxids_l) == 1:
                    #     sscinames = sscinames_l[0]
                    #     qtaxid = staxids_l[0]
                    # # if not the first item and id different from before: get name
                    # elif i != 0 and id_before != id_now:
                    #     sscinames = self.ids.ncbi_parser.get_name_from_id(qtaxid)
                    # elif i == 0:  # for first item in redundant data, always get info
                    #     sscinames = self.ids.ncbi_parser.get_name_from_id(qtaxid)
                    # else:  # if id_before and id_now were the same, we do not need to add same seq again
                    #     continue
                    try:
                        # assert str(staxids) in staxids_l, (str(staxids), staxids_l)
                        # next vars are used to stop loop if all taxids were found
                        found_taxids.add(staxids)
                        found_spn.add(sscinames)
                        if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                            query_dict[gb_acc] = \
                                {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                 'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                 'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                    except AssertionError:
                        sys.stderr.write("Taxon Id {} was not initally found in the blast search - corresponding Genbank accession number is: {}. Sequence will not be added. Check with ncbi. \n".format(staxids, gb_acc))
         
            # same loop as above, only that it does more blastdbcmd's
            # this is used as sometimes the if above does not yield in stop_while == True,
            # through different taxa names
            elif count >= 1 and stop_while is False:
                # debug("count>1")
                for i in range(0, len(sallseqid_l)):
                    if len(found_taxids) == len(staxids_l):
                        break
                    gb_acc = get_acc_from_blast(sallseqid_l[i])
                    gi_id = get_gi_from_blast(sallseqid_l[i])
                    stitle = salltitles_l[i]
                    # if gb acc was already read in before stop the for loop
                    if gb_acc in query_dict or gb_acc in self.data.gb_dict:
                        stop_while = True
                        break
                 
                    staxids = self.get_taxid_from_acc(gb_acc)[0]
                    sscinames = self.ids.ncbi_parser.get_name_from_id(staxids)
                    assert str(staxids) in staxids_l, (staxids, staxids_l)
                    # next vars are used to stop loop if all taxids were found
                    found_taxids.add(staxids)
                    found_spn.add(sscinames)
                    if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                        query_dict[gb_acc] = \
                             {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                              'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                              'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
        return query_dict

    def read_local_blast_query(self, fn_path):
        """ Implementation to read in results of local blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gb_dict dictionaries
        """
        debug("read_local_blast_query")
        query_dict = {}
        with open(fn_path, mode="r") as infile:
            for lin in infile:
                sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, salltitles, sallseqid = lin.strip().split('\t')
                gb_acc = get_acc_from_blast(sseqid)
                gi_id = get_gi_from_blast(sseqid)                
                sseq = sseq.replace("-", "") #TODO here is where we want to grab the full sequence MK: I wrote a batch query for the seqs we are interested. Makes it faster.
                sscinames = sscinames.replace(" ", "_").replace("/", "_")
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                stitle = salltitles
                # get additional info only for seq that pass the eval
                if evalue < float(self.config.e_value_thresh):
                    # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
                    # Get all of them! (they can be of a different taxon ids = get redundant seq info)
                    if len(sallseqid.split(";")) > 1:
                        var_list = [gb_acc, gi_id, sseq, staxids, sscinames, pident, evalue, bitscore, stitle, sallseqid]
                        query_dict = self.get_new_seqs_for_mergedseq(var_list, query_dict)
                        taxids_l = staxids.split(";")
                    else:  # if there are no non-redundant data
                        staxids = int(staxids)
                        self.ids.spn_to_ncbiid[sscinames] = staxids
                        if gb_acc not in self.ids.acc_ncbi_dict:  # fill up dict with more information.
                            self.ids.acc_ncbi_dict[gb_acc] = staxids
                        if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                            query_dict[gb_acc] = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                                  'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                                  'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
        # debug("key in query")
        for key in query_dict.keys():
            if float(query_dict[key]["evalue"]) < float(self.config.e_value_thresh):
                gb_acc = query_dict[key]["accession"]
                if len(gb_acc.split(".")) >= 2:
                    # skip ones we already wanted to add earlier, exception if add_lower_taxa
                    if gb_acc not in self.data.gb_dict or self.config.add_lower_taxa is True:
                        # debug("add gb_dict")
                        self.new_seqs[gb_acc] = query_dict[key]["sseq"]
                        self.data.gb_dict[gb_acc] = query_dict[key]
                # else:
                    # debug("was added before")
            else:
                fn = open("{}/blast_threshold_not_passed.csv".format(self.workdir), "a+")
                fn.write("blast_threshold_not_passed: {}, {}, {}\n".format(query_dict[key]["sscinames"], query_dict[key]["accession"],
                         query_dict[key]["evalue"]))
                fn.close()

    def get_taxid_from_acc(self, gb_acc):
        """
        Use the blastdb to get the taxon id from a queried gb acc. 

        Sometimes there are more than a singel id, as database has merged redundant seqs.

        :param gb_acc: Genbank accession number
        :return: list of taxon ids associated with the GB id - there are multiple because of the merging of redundant data
        """
        if not os.path.exists("{}/tmp".format(self.workdir)):
            os.mkdir("{}/tmp".format(self.workdir))
        fn = "{}/tmp/tmp_search.csv".format(self.workdir)
        fn_open = open(fn, "w+")                                        
        fn_open.write("{}\n".format(gb_acc))
        fn_open.close()
        cmd1 = "blastdbcmd -db {}/nt  -entry_batch {} -outfmt %T -out {}/tmp/tax_id_{}.csv".format(self.config.blastdb, fn, self.workdir, gb_acc)
        if not self.config.blastdb == "./tests/data/precooked/testing_localdb":
            os.system(cmd1)
            f = open("{}/tmp/tax_id_{}.csv".format(self.workdir, gb_acc))
            tax_id_l = []
            for line in iter(f):
                line = line.rstrip().lstrip()
                tax_id_l.append(int(line))
        else: 
            f = open("./tests/data/precooked/testing_localdb/tax_id_{}.csv".format(gb_acc))
            tax_id_l = []
            for line in iter(f):
                line = line.rstrip().lstrip()
                tax_id_l.append(int(line))
        f.close() 
        return tax_id_l

    def read_unpublished_blast_query(self):
        """
        Reads in the blast files generated during local_blast_for_unpublished() and adds seq to self.data.gb_dict and
        self.new_seqs.

        """
        debug("read unpublished blast query")
        output_blast = "output_tst_fn.xml"
        gb_counter = 1
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        # with cd(os.path.join(self.workdir, "blast")):
        xml_file = open(output_blast)
        os.chdir(general_wd)
        blast_out = NCBIXML.parse(xml_file)
        fn = open("{}/not_added_local_seq.csv".format(self.workdir), "a")
        fn.write("not_added_local_seq")
        for blast_record in blast_out:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    local_id = alignment.title.split("|")[-1].split(" ")[-1]
                    if float(hsp.expect) < float(self.config.e_value_thresh):
                        if local_id not in self.data.gb_dict:  # skip ones we already have
                            unpbl_local_id = "unpubl_{}".format(local_id)
                            self.new_seqs[unpbl_local_id] = hsp.sbjct
                            # debug(self.new_seqs[unpbl_local_id])
                            self.data.gb_dict[unpbl_local_id] = {'title': "unpublished", 'localID': local_id}
                            self.data.gb_dict[unpbl_local_id].update(
                                self.data.unpubl_otu_json['otu{}'.format(local_id.replace("_", "").replace("-", ""))])
                            gb_counter += 1
                    else:
                        fn.write("{}: {}".format(alignment.title.split("|")[-1].split(" ")[-1], hsp.expect))
                        writeinfofiles.write_not_added_info(self, local_id, "threshold not passed")
                        # needs to be deleted from gb_dict,
                        # maybe we find a better fitting blast query seq and then it might get added
                        del self.data.gb_dict[unpbl_local_id]  # print(some)
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from unpublished database\n".format(len(self.new_seqs)))

    def read_webbased_blast_query(self, fn_path):
        """ Implementation to read in results of web blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gb_dict dictionaries
        """
        result_handle = open(fn_path)
        try:
            if _VERBOSE:
                sys.stdout.write(".")
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if float(hsp.expect) < float(self.config.e_value_thresh):
                            gb_id = alignment.title.split("|")[3]  # 1 is for gi
                            if len(gb_id.split(".")) == 1:
                                debug(gb_id)
                            if gb_id not in self.data.gb_dict:  # skip ones we already have
                                # gb_id = int(alignment.title.split('|')[1])  # 1 is for gi
                                # assert type(gb_id) is int
                                # SHOULD NOT BE NECESSARY....IS WEBBLAST HAS THE TAXON ALREADY LIMITED
                                # if len(self.acc_list_mrca) >= 1 and (gb_id not in self.acc_list_mrca):
                                #     pass
                                # else:
                                taxid,taxname, seq = self.ids.get_tax_seq_acc(gb_id)
                                self.new_seqs[gb_id] = seq
                                gi_id = alignment.title.split('|')[1]
                                gb_acc = alignment.__dict__['accession']
                                stitle = alignment.__dict__['title']
                                hsps = alignment.__dict__['hsps']
                                length = alignment.__dict__['length']
                                query_dict = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'title': stitle,
                                              'length': length, 'hsps': hsps}
                                self.data.gb_dict[gb_id] = query_dict
                        else:
                            writeinfofiles.write_not_added_info(self, gb_acc, "threshold not passed")
                            # needs to be deleted from gb_dict,
                            # maybe we find a better fitting blast query seq and then it might get added
                            del self.data.gb_dict[gb_acc]
        except ValueError:
            sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

    def read_blast_wrapper(self, blast_dir=None):
        """reads in and processes the blast xml files

        :param blast_dir: path to directory which contains blast files
        :return: fills different dictionaries with information from blast files
        """
        debug("read_blast_wrapper")
        if blast_dir:
            if _VERBOSE:
                sys.stdout.write("blast dir is {}\n".format(blast_dir))
            self.blast_subdir = os.path.abspath(blast_dir)
        else:
            if _VERBOSE:
                sys.stdout.write("blast dir is {}\n".format(self.blast_subdir))
            if not os.path.exists(self.blast_subdir):
                os.mkdir(self.blast_subdir)
        if self.unpublished:
            self.read_unpublished_blast_query()
        else:
            if not self._blasted:
                self.run_blast_wrapper()
            assert os.path.exists(self.blast_subdir)
            for taxon in self.data.aln:
                # debug(self.config.blast_loc)
                fn = None
                if self.config.blast_loc == "local":
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                if self.config.gb_id_filename is True: #TODO what is this doing?
                    fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label) 
                    if fn is None:
                        fn = self.data.otu_dict[taxon.label].get('^user:TaxonName', taxon.label)
                    fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                else:
                    fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                if _DEBUG:
                    sys.stdout.write("attempting to read {}\n".format(fn_path))
                if os.path.isfile(fn_path):
                    if self.config.blast_loc == 'local':  # new method to read in txt format
                        self.read_local_blast_query(fn_path)
                    else:
                        self.read_webbased_blast_query(fn_path)
        self.date = str(datetime.date.today())
#        debug("len new seqs dict after evalue filter")
#        debug(len(self.new_seqs))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from GenBank after evalue filtering\n".format(len(self.new_seqs)))

        self._blast_read = 1

    def get_sp_id_of_otulabel(self, label): 
        #TODO This was doing a bunch of searches that sould have already happened when creating the OTU.
        #If we didn't find an ncbi_id then, we shouldn't look again.

        """Get the species name and the corresponding ncbi id of the otu.

        :param label: otu_label = key from otu_dict
        :return: ncbi id of corresponding label
        """
        # debug("get_tax_id_of_otulabel")
        ncbi_id = self.data.otu_dict[label].get("^ncbi:taxon",0)
#        debug("otu {} has taxon id {}".format(label, ncbi_id))
        return ncbi_id


    def seq_dict_build(self, seq, label, seq_dict):
        """takes a sequence, a label (the otu_id) and a dictionary and adds the
        sequence to the dict only if it is not a subsequence of a
        sequence already in the dict.
        If the new sequence is a super sequence of one in the dict, it
        removes that sequence and replaces it

        :param seq: sequence as string, which shall be compared to existing sequences
        :param label: otu_label of corresponding seq
        :param seq_dict: the tmp_dict generated in add_otu()
        :return: updated seq_dict
        """
        if label == None: #in case of add_otu failure, doean't edit dict 
            debug("otu_id None")
            sys.stderr.write("otu_id None")
            return seq_dict
        id_of_label = self.get_sp_id_of_otulabel(label)
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        continue_search = False
        never_add = False
        for tax_lab in tax_list:
            existing_id = self.get_sp_id_of_otulabel(tax_lab)
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(new_seq) >= sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * self.config.maxlen:
                debug("seq not added because it's too long...")
                self.data.otu_dict[label]['^physcraper:status'] = "not added; sequence is too long"
                gb_id = self.data.otu_dict[label]["^ncbi:accession"]
                # if gb_id not in self.gb_not_added:
                #     self.gb_not_added.append(gb_id)
                #TODO merge with other label search
                if id_of_label in self.ids.ncbiid_to_spn.keys():
                    tax_name = self.ids.ncbiid_to_spn[id_of_label]
                else:
                    tax_name = self.ids.ncbi_parser.get_name_from_id(id_of_label)
                cutoff = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * self.config.maxlen
                reason = "sequence too long: new seq len ({}) vs.  cutoff ({})".format(len(new_seq), cutoff)
                writeinfofiles.write_not_added(id_of_label, tax_name, gb_id, reason, self.workdir)
                # writeinfofiles.write_not_added_info(self, local_id, "threshold not passed")
                # needs to be deleted from gb_dict,
                # maybe we find a better fitting blast query seq and then it might get added
#                del self.data.gb_dict[gb_id]
                                    
            elif len(inc_seq) >= len(new_seq):  # if seq is identical and shorter
                if inc_seq.find(new_seq) != -1:
                    if type(existing_id) == int and existing_id != int(id_of_label):  # different otus, add
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, "
                                             "but different species name\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; " \
                                                                          "subsequence, but different taxon"
                        seq_dict[label] = seq
                        debug("{} and {} are subsequences, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:  # subseq of same otu
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        debug("{} not added, subseq of {}".format(id_of_label, existing_id))
                        never_add = True
                        continue
            else:  # if seq is longer and identical
                if new_seq.find(inc_seq) != -1:
                    if self.data.otu_dict[tax_lab].get('^physcraper:status') == "original":
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of original seq {}, "
                                             "both kept in alignment\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        debug("{} and  {} added".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    elif type(existing_id) == int and existing_id != id_of_label:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, but different "
                                             "taxon\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; supersequence, " \
                                                                          "but different taxon"
                        seq_dict[label] = seq
                        debug("{} and  {} supersequence, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    elif self.data.otu_dict[tax_lab].get('^physcraper:status') == "local seq":
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of local seq {}, "
                                             "both kept in alignment\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        debug("{} and  {} added".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, {} added "
                                             "and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added in place of {}".format(tax_lab)
                        self.data.otu_dict[tax_lab]['^physcraper:status'] = "deleted, {} is supersequence".format(label)
                        debug("{} added, instead of  {}".format(id_of_label, existing_id))
                        self.del_superseq.add(tax_lab)
                        continue_search = True
                        continue

        if continue_search is True or never_add is True:
            if (self.data.otu_dict[label]['^physcraper:status'].split(' ')[0] in self.seq_filter) or never_add is True:
                if label in seq_dict.keys():
                    del seq_dict[label]
                if label in self.data.aln.taxon_namespace or label in self.data.tre.taxon_namespace:
                    self.data.remove_taxa_aln_tre(label)
                # should not be necessary, information what happened to seq should have been added in lines before
                # else:
                #     debug("label was never added to aln or tre")
                # Note: should not be the word 'deleted', as this is used in self.seq_filter
                # self.data.otu_dict[label]['^physcraper:status'] = "removed in seq dict build"
                return seq_dict
        if _VERBOSE:
            sys.stdout.write(".")
            if i % 50 == 0:
                sys.stdout.write("\n")
        seq_dict[label] = seq
        return seq_dict


    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        """
        debug("remove identical seqs")
        if len(self.new_seqs_otu_id) > 0:
            if _DEBUG:
                sys.stdout.write("running remove identical twice in a row"
                                 "without generating new alignment will cause errors. skipping\n")
            return
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        # Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        # need to re-calculate orig_seq_len before using it
        self.data.orig_seqlen = [len(self.data.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in
                                 self.data.aln]
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        assert self.config.seq_len_perc <= 1, \
            ("your config seq_len_param is not smaller than 1: {}".format(self.config.seq_len_perc))
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        self.del_superseq = set()  # will contain deleted superseqs for the assert below
        # all_added_gi is to limit the adding to new seqs, if we change the rank of filtering later
        all_added_gi = set()
        if self.config.add_lower_taxa is True:
            for key in self.data.otu_dict.keys():
                # debug(self.data.otu_dict[key]['^physcraper:status'].split(' ')[0])
                if self.data.otu_dict[key]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    # debug(self.data.otu_dict[key].keys())
                    if "^ncbi:accession" in self.data.otu_dict[key]:
                        all_added_gi.add(self.data.otu_dict[key]["^ncbi:accession"])
        for gb_id, seq in self.new_seqs.items():
            added = False
            assert gb_id in self.data.gb_dict.keys(), (gb_id, self.data.gb_dict.keys())
            if gb_id not in all_added_gi:
                reason = ""
                all_added_gi.add(gb_id)
                if len(gb_id.split(".")) == 1:
                    debug("problem gb_id {}".format(gb_id))
                if self.blacklist is not None and gb_id in self.blacklist:
                    debug("gb_id in blacklist, not added")
                    pass
                elif gb_id in self.newseqs_acc:  # added to increase speed. often seq was found in another blast file
                    debug("passed, was already added")
                    pass
                else:
                    # debug(len(seq.replace("-", "").replace("N", "")))
                    # debug(seq_len_cutoff)
                    if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                        # debug(self.config.blast_loc)
                        if self.config.blast_loc != "remote":
                            # debug("ADD")
                            tax_name = None
                            # ######################################################
                            # ### new implementation of rank for delimitation
                            # if type(self.mrca_ncbi) is int:
                            #     mrca_ncbi = self.mrca_ncbi
                            #     debug("int")
                            #     debug(mrca_ncbi)
                            # elif len(self.mrca_ncbi) == 1:
                            #     mrca_ncbi = list(self.mrca_ncbi)[0]
                            #     debug(mrca_ncbi)
                            # else:
                            #     debug(self.mrca_ncbi)
                            #     debug("should not happen. think about something to do!")
                            #     mrca_ncbi = None
                            # # debug(stop)
                            # # debug(mrca_ncbi)
                            # #############################################

                            # rank_mrca_ncbi = self.ids.ncbi_parser.get_rank(mrca_ncbi)
                            # get rank to delimit seq to ingroup_mrca
                            ncbi_id, tax_name = get_tax_info_from_acc(gb_id, self.data, self.ids)

                            rank_mrca_ncbi = set()
                            for mrca in self.mrca_ncbi_list:
                                #debug("mrca")
                                #debug(mrca)
                                rank_mrca = self.ids.ncbi_parser.get_rank(mrca)
                                rank_mrca_ncbi.add(rank_mrca)
                                try:  # sometimes ncbi has wrong id linked
                                    input_rank_id = self.ids.ncbi_parser.match_id_to_mrca(ncbi_id, mrca) 
                                    #ToDO changing this funtion to return T if is is within mrca and F if it is not will imporve readability
                                    #(I don't think we use the actual value agin)
                                except IndexError:  # this is for the wrong ncbi link, get right tax_id and do search again
                                    debug("wrong tax_id given by ncbi? tax_id: {}, tax_name:{}".format(ncbi_id, tax_name))
                                    old_ncbi_id = ncbi_id
                                    ncbi_id = self.ids.ncbi_parser.get_id_from_name(tax_name)
                                    # debug(ncbi_id)
                                    self.data.gb_dict[gb_id]['staxids'] = ncbi_id
                                    input_rank_id = self.ids.ncbi_parser.match_id_to_mrca(ncbi_id, mrca)
                                    sys.stderr.write("accession number {} was assigned taxon id {}, which is not in ncbi taxonomy\n".format(gb_id, old_ncbi_id))
                                    sys.stderr.write("searched on name {} and assigned taxon id {} instead \n".format(tax_name, ncbi_id))
                                #debug([input_rank_id, mrca])
                                if int(input_rank_id) == int(mrca):  # belongs to ingroup mrca -> add to data
                                   # debug("input belongs to same mrca")
                                    self.newseqs_acc.append(gb_id)
                                    otu_id = self.data.add_otu(gb_id, self.ids)
                                    self.seq_dict_build(seq, otu_id, tmp_dict)
                                    added = True
                                    # debug(some)
                                else:
                                # elif gb_id not in self.gb_not_added:  # if not, leave it out
                                    debug("does not match")
                                    # debug([mrca, input_rank_id])
                                    reason = "not_part_of_mrca: {}".format(self.mrca_ncbi_list)
                        else:
                            added = True
                            self.newseqs_acc.append(gb_id)
                            otu_id = self.data.add_otu(gb_id, self.ids)
                            self.seq_dict_build(seq, otu_id, tmp_dict)
                    else:
                        len_seq = len(seq.replace("-", "").replace("N", ""))
                        reason = "seqlen_threshold_not_passed: len seq: {}, min len: {}".format(len_seq, seq_len_cutoff)
 
                # sequences that could not be added are written out to file
                if added is False:
                    if "sscinames" in self.data.gb_dict[gb_id]:
                            tax_name = self.data.gb_dict[gb_id]['sscinames']
                            ncbi_id = self.data.gb_dict[gb_id]['staxids']
                    else:
                        tax_name = None
                        ncbi_id = None
                    writeinfofiles.write_not_added(ncbi_id, tax_name, gb_id, reason, self.workdir)

        # this assert got more complicated, as sometimes superseqs are already deleted in seq_dict_build() ->
        # then subset assert it not True
        old_seqs_ids = set()
        for tax in old_seqs:
            old_seqs_ids.add(tax)
        tmp_dict_plus_super = set()
        for item in tmp_dict.keys():
            tmp_dict_plus_super.add(item)
        for item in self.del_superseq:
            tmp_dict_plus_super.add(item)
        # assert old_seqs_ids.issubset(tmp_dict.keys()), ([x for x in old_seqs_ids if x not in tmp_dict.keys()])
        assert old_seqs_ids.issubset(tmp_dict_plus_super), ([x for x in old_seqs_ids if x not in tmp_dict_plus_super])
        for tax in old_seqs:
            if tax in tmp_dict:
                del tmp_dict[tax]
        self.new_seqs_otu_id = tmp_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        debug("len new seqs dict after remove identical")
        debug(len(self.new_seqs_otu_id))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from Genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def dump(self, filename=None, recursion=100000):
        """writes out class to pickle file.
        We need to increase the recursion depth here, as it currently fails with the standard run.

        :param filename: optional filename
        :param recursion: pickle often failed with recursion depth, that's why it's increased here
        :return: writes out file
        """
        current = sys.getrecursionlimit()
        sys.setrecursionlimit(recursion)

        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/scrape_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi, pickle.HIGHEST_PROTOCOL)
        sys.setrecursionlimit(current)

    def write_query_seqs(self):
        """writes out the query sequence file"""
        debug("write query seq")
        if not self._blast_read:
            self.read_blast_wrapper()
        self.newseqs_file = "{}.fasta".format(self.date)
        fi = open("{}/{}".format(self.workdir, self.newseqs_file), "w")
        if _VERBOSE:
            sys.stdout.write("writing out sequences\n")
        for otu_id in self.new_seqs_otu_id.keys():
            if otu_id not in self.data.aln:  # new seqs only
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1

    def align_query_seqs(self, papara_runname="extended"):
        """runs papara on the tree, the alignment and the new query sequences

        :param papara_runname: possible file extension name for papara
        :return: writes out files after papara run/aligning seqs
        """
        cwd = os.getcwd()
        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        if _VERBOSE:
            sys.stdout.write("aligning query sequences \n")
        self.data._reconcile()  # I think reconcile is what was needed here...instead of alien hack
        # note: sometimes there are still sp in any of the aln/tre
        # hack for the alien taxa thing
        self.remove_alien_aln_tre()
        self.data.write_papara_files()
        os.chdir(self.workdir)  # Clean up dir moving
        # with cd(self.workdir):
        assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
        try:
            subprocess.check_call(["papara",
                                   "-t", "random_resolve.tre",
                                   "-s", "aln_ott.phy",
                                   #  "-j", "{}".format(self.config.num_threads),  # FIXME: Does not work on some machines
                                   "-q", self.newseqs_file,
                                   "-n", papara_runname])  # FIXME directory ugliness
            if _VERBOSE:
                sys.stdout.write("Papara done")
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, grepexc.output)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it installed?\n")
                sys.exit(-5)
            # handle file not found error.
            else:
                # Something else went wrong while trying to run `wget`
                raise
        path = "{}/papara_alignment.{}".format(self.workdir, papara_runname)
        assert os.path.exists(path), "{path} does not exists".format(path=path)
        os.chdir(cwd)
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment."
                                                    "{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        self._query_seqs_aligned = 1

    def remove_alien_aln_tre(self):
        """Sometimes there were alien entries in self.tre and self.aln.

        This function ensures they are properly removed."""

        treed_tax = set()
        for leaf in self.data.tre.leaf_nodes():
            treed_tax.add(leaf.taxon)
        aln_tax = set()
        for tax, seq in self.data.aln.items():
            aln_tax.add(tax)
        prune = treed_tax ^ aln_tax
        del_tre = []
        del_aln = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_tax)
            if taxon in aln_tax:
                # debug(taxon)
                del_aln.append(taxon)
            if taxon in treed_tax:
                del_tre.append(taxon)
        # debug(del_aln)
        # debug(del_tre)
        self.data.aln.remove_sequences(del_aln)
        self.data.tre.prune_taxa(del_tre)
        for tax_lab in self.data.aln.taxon_namespace:
            if tax_lab not in self.data.tre.taxon_namespace:
                sys.stderr.write("tax {} not in tre. This is an alien name in the data.\n".format(tax_lab))
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_lab in self.data.tre.taxon_namespace:
            if tax_lab not in self.data.aln.taxon_namespace:
                sys.stderr.write("tax {} not in aln. This is an alien name in the data.\n".format(tax_lab))
                self.data.remove_taxa_aln_tre(tax_lab)
        # this should not need to happen here
        # self.data.prune_short()
        # self.data.trim()

    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree."""
        if self.backbone is True:
            with cd(self.workdir):
                backbonetre = Tree.get(path="{}/backbone.tre".format(self.workdir),
                                       schema="newick",
                                       preserve_underscores=True)

                backbonetre.resolve_polytomies()
                backbonetre.write(path="random_resolve.tre", schema="newick", unquoted_underscores=True)

        if os.path.exists("RAxML_labelledTree.PLACE"):
            os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
        if _VERBOSE:
            sys.stdout.write("placing query sequences \n")
        with cd(self.workdir):
            try:
                debug("try")
                subprocess.call(["raxmlHPC-PTHREADS",
                                 "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-f", "v",
                                 "-s", "papara_alignment.extended",
                                 "-t", "random_resolve.tre",
                                 "-n", "PLACE"])
                placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                    schema="newick",
                                    preserve_underscores=True)
            except:
                try:
                    subprocess.call(["raxmlHPC",
                                     "-m", "GTRCAT",
                                     "-f", "v",
                                     "-s", "papara_alignment.extended",
                                     "-t", "random_resolve.tre",
                                     "-n", "PLACE"])
                    placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                        schema="newick",
                                        preserve_underscores=True)
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        sys.stderr.write("failed running raxmlHPC. Is it installed?")
                        sys.exit(-6)
                    # handle file not
                    # handle file not found error.
                    else:
                        # Something else went wrong while trying to run `wget`
                        raise
            placetre.resolve_polytomies()
            for taxon in placetre.taxon_namespace:
                if taxon.label.startswith("QUERY"):
                    taxon.label = taxon.label.replace("QUERY___", "")
            placetre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True)
        self._query_seqs_placed = 1

    def est_full_tree(self, path="."):
        """Full raxml run from the placement tree as starting tree.
        The PTHREAD version is the faster one, hopefully people install it. if not it falls back to the normal raxml.
        the backbone options allows to fix the sceleton of the starting tree and just newly estimates the other parts.
        """
        debug("est full tree")
        cwd = os.getcwd()
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        try:
            num_threads = int(self.config.num_threads)
            if self.backbone is not True:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
            else:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")

            if self.backbone is not True:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
            else:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
        os.chdir(cwd)
        self._full_tree_est = 1

    def calculate_bootstrap(self):
        """Calculates bootstrap and consensus trees.

        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        -z: specifies file with multiple trees

        """
        debug("calculate bootstrap")
        cwd = os.getcwd()
        os.chdir(self.workdir)
        # with cd(self.workdir):
        # # check if job was started with mpi
        # # this checks if the whole file was started as mpiexec
        # env_var = [os.environ.get('PMI_RANK'), os.environ.get('PMI_SIZE'), os.environ.get('OMPI_COMM_WORLD_SIZE')]
        # mpi = False
        # for var in env_var:
        #     if var is not None:
        #         mpi = True
        # check if job was started with mpi
        # this checks if actual several cores and nodes were allocated
        ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
        nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
        # env_var = int(nnodes) * int(ntasks)
        mpi = False
        if nnodes is not None and ntasks is not None:
            env_var = int(nnodes) * int(ntasks)
            mpi = True
        if mpi:
            debug("run with mpi")
            subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2",
                             # "raxmlHPC-PTHREADS", "-T", "{}".format(num_threads),
                             "-m", "GTRCAT",
                             "-s", "previous_run/papara_alignment.extended",
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "{}".format(self.date)])
        else:
            try:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-s", "previous_run/papara_alignment.extended",
                                 "-p", "1", "-b", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])
            except:
                subprocess.call(["raxmlHPC",
                                 "-m", "GTRCAT",
                                 "-s", "previous_run/papara_alignment.extended",
                                 "-p", "1", "-b", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])

        try:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-s", "previous_run/papara_alignment.extended",
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "all{}".format(self.date)])

            # strict consensus:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "STRICT",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "StrictCon{}".format(self.date)])
            # majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MR",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "MR_{}".format(self.date)])
            # extended majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MRE",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "EMR{}".format(self.date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")
            # make bipartition tree
            # is the -f b command
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-s", "previous_run/papara_alignment.extended",
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "all{}".format(self.date)])
            # strict consensus:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "STRICT",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "StrictCon{}".format(self.date)])
            # majority rule:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "MR",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "MR_{}".format(self.date)])
            # extended majority rule:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "MRE",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "EMR{}".format(self.date)])
        os.chdir(cwd)

    def remove_blacklistitem(self):
        """This removes items from aln, and tree, if the corresponding Genbank identifer were added to the blacklist.

        Note, that seq that were not added because they were similar to the one being removed here, are lost
        (that should not be a major issue though, as in a new blast_run, new seqs from the taxon can be added.)
        """
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            acc = self.data.otu_dict[tax.label].get("^ncbi:accession")
            if gi_id in self.blacklist or acc in self.blacklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, Genbank identifier is part of blacklist"
        # this should not need to happen here: prune_short; instead...
        self.data.check_tre_in_aln()
        # self.data.prune_short()
        # debug(self.data.tre.as_string(schema='newick'))

    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignment with the expanded ones"""
        debug("generate streamed aln")
        if self.blacklist:
            self.remove_blacklistitem()
        debug(len(self.new_seqs))
        debug(len(self.new_seqs_otu_id))
        if len(self.new_seqs) == 0 or len(self.new_seqs_otu_id) == 0:
            if _VERBOSE:
                sys.stdout.write("No new sequences found.\n")
            # self.repeat = 0
            self.calculate_final_tree()
            self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))
        elif len(self.new_seqs) > 0:
            self.data.write_files()  # should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:  # TODO rename to something more intuitive
                self.data.check_tre_in_aln()
                self.write_query_seqs()
                self.align_query_seqs()
                self.place_query_seqs()
                self.data.prune_short()
                self.data.trim()
                self.est_full_tree()
                self.data.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.workdir, self.date),
                                         schema="newick",
                                         preserve_underscores=True,
                                         taxon_namespace=self.data.aln.taxon_namespace)
                self.data.write_files()
                if os.path.exists("{}/previous_run".format(self.workdir)):
                    prev_dir = "{}/previous_run{}".format(self.workdir, self.date)
                    i = 0
                    while os.path.exists(prev_dir):
                        i += 1
                        prev_dir = "{}/previous_run{}".format(self.workdir, self.date) + str(i)
                    os.rename("{}/previous_run".format(self.workdir), prev_dir)
                if self.config.gb_id_filename is not True:
                    os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    if not os.path.exists("{}/previous_run".format(self.workdir)):
                        os.makedirs('{}/previous_run/'.format(self.workdir))
                    if self.config.gb_id_filename is not True:
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                    else:
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file),
                          "{}/previous_run/newseqs.fasta".format(self.workdir))
                self.data.write_labelled(label='^physcraper:TaxonName', add_gb_id=True)
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {}  # Wipe for next run
                self.new_seqs_otu_id = {}
                self.repeat = 1
            else:
                self.calculate_final_tree()
                self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))
        #     else:
        #         if _VERBOSE:
        #             sys.stdout.write("No new sequences after filtering.\n")
        #         # self.repeat = 0
        #         self.calculate_final_tree()
        #         self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))

        # else:
        #     if _VERBOSE:
        #         sys.stdout.write("No new sequences found.\n")
        #     # self.repeat = 0
        #     self.calculate_final_tree()
        #     self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))

        self.reset_markers()

        filter_by_local_blast.del_blastfiles(self.workdir)  # delete local blast db
        self.data.dump()
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))

    def calculate_final_tree(self):
        """Calculates the final tree using a trimmed alignment.

        :return: final PS data
        """
        debug("calculate final tree")
        self.data.write_files(treepath="physcraper_final_notrim.tre", alnpath="physcraper_final_notrim.fas")
        self.data.prune_short()
        self.data.trim()
        self.data.write_files(treepath="physcraper_final_trim.tre", alnpath="physcraper_final_trim.fas")
        if os.path.exists("[]/previous_run".format(self.workdir)):
            self.est_full_tree(path="previous_run")
        else:
            self.est_full_tree()
        self.repeat = 0
        self.calculate_bootstrap()

    def write_unpubl_blastdb(self, path_to_local_seq):
        """Adds local sequences into a  local blast database, which then can be used to blast aln seq against it
        and adds sequences that were found to be similar to input.
        If this option is used, it queries against local database first and only in "2" round
        it goes back to blasting against GenBank

        :param path_to_local_seq: path to the local seqs that shall be added
        :return: writes local blast databases for the local sequences
        """
        debug("add_local_seq")
        self.path_to_local_seq = path_to_local_seq
        localfiles = os.listdir(path_to_local_seq)
        for index, item in enumerate(localfiles):
            item = str(item)
            if item.startswith(".~"):
                localfiles[index] = None
        localfiles = filter(None, localfiles)
        for filename in localfiles:
            filepath = "{}/{}".format(path_to_local_seq, filename)
            open_file = open(filepath)
            content = open_file.readlines()
            content = [x.strip() for x in content]
            content = filter(None, content)  # fastest
            count = 0
            gb_id_l = content[::2]
            seq_l = content[1::2]
            # in current setup 1 seq per file, but this is written in a way,
            # that a file with multiple seqs can be read in as well
            for i in xrange(0, len(gb_id_l)):
                key = gb_id_l[i].replace(">", "")
                count = count + 1
                seq = seq_l[i]
                filter_by_local_blast.write_filterblast_db(self.workdir, key, seq, fn="local_unpubl_seq")
        with cd(os.path.join(self.workdir, "blast")):
            cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
            os.system(cmd1)


class FilterBlast(PhyscraperScrape):
    """Takes the Physcraper Superclass and filters the ncbi blast results to only include a subset of the sequences.

    They can be filtered by number or by rank and number. The feature can be useful for non-population-level studies,
    e.g. analyses which require a single representative per taxon (threshold = 1) or to check the monophyly of taxa
    without having to deal with over-representation of few taxa (e.g. threshold = 4, which allows to get a good overview
    of what is available without having some taxa being represented by high numbers of sequences).
    The second option (downtorank), which is optional, allows to filter according to taxonomic levels, e.g. getting
    a number of representative sequences for a genus or lineage it can also be used to deal with subspecies.

    Existing self objects are:

        self.sp_d: dictionary

                key = species name/id
                value = list with otu_dict entry

        self.sp_seq_d: dictionary

                key = species name/id
                value = dictionary (Is overwritten every 'round')

                    key = otuID
                    value = seq.
        self.filtered_seq: dictionary. Is used as the self.new_seqs equivalent from Physcraper, just with fewer seqs. Is overwritten every 'round'

                key = otuID,
                val = seq.
        self.downtorank: optional string defining the level of taxonomic filtering, e.g. "species", "genus"
    """

    # TODO MK: self.sp_d = {} does not need to contain all otu_dict info, key is sufficient

    def __init__(self, data_obj, ids, ingroup_mrca=None):
        super(FilterBlast, self).__init__(data_obj, ids, ingroup_mrca)
        debug("start derived class init")
        # additional things that are needed for the filtering process
        debug(self.mrca_ncbi_list)
        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}
        self.downtorank = None
        self.threshold = None

    def add_setting_to_self(self, downtorank, threshold):
        """
        Add FilterBlast items to self.

        Currently used by some wrapper functions.

        :param downtorank: rank which defines your level of Filtering
        :param threshold: number, defining how many seq per rank do you want to keep
        :return:
        """
        self.threshold = threshold
        self.downtorank = downtorank

    def sp_dict(self, downtorank=None):
        """Takes the information from the Physcraper otu_dict and makes a dict with species id as key and
        the corresponding seq information from aln and blast seq, it returns self.sp_d.

        This is generated to make information for the filtering class more easily available. self.sp_d sums up which
        information are available per taxonomic concept and which have not already been removed during either
        the remove_identical_seq steps or during a filtering run of an earlier cycle.

        Note: has test, test_sp_d.py

        :param downtorank: string defining the level of taxonomic filtering, e.g. "species", "genus"
        :return: self.sp_d
        """
        #Edited to pass in list of otu_ids rather than set of dictionaries, to make getting squence by id easier in sp_seq_d
        self.downtorank = downtorank
        debug("make sp_dict")
        self.sp_d = {}
        for otu_id in self.data.otu_dict:
            if self.data.otu_dict[otu_id]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                #TODO these lines should not be necessary, as we have already mapped the ncbi_id when creating the otu
                #can replace with 
                #tax_id = self.get_sp_id_of_otulabel(otu_id)
                tax_name = self.ids.find_name_otu(otu_dict_entry=self.data.otu_dict[otu_id])
                if len(tax_name.split("(")) > 1:
                    tax_name = tax_name.split("(")[0]
                tax_name = str(tax_name).replace(" ", "_")
                if self.config.blast_loc == 'remote':
                    sys.stderr.write("Filtering by taxon not functional for remote ncbi searches yet.")
                    sys.exit(-7)
                    #if '^ncbi:accession' in self.data.otu_dict[otu_id]:
                    #     gb_id = self.data.otu_dict[otu_id]['^ncbi:accession']
                    #    if len(gb_id.split(".")) == 1:
                    #        debug(gb_id)
                    #    if gb_id in self.ids.acc_ncbi_dict:
                    #        tax_id = self.ids.acc_ncbi_dict[gb_id]
                    #tax_id = self.ids.get_rank_info_from_web(taxon_name=tax_name)
                    # print(tax_name)
                    # print(self.ids.otu_rank.keys())
                    # tax_id = self.ids.otu_rank[tax_name]["taxon id"]
                else:
                    #We already earched for the taxon_id when creating the otu dict entry.
                    # We want to make sure that the ncbi_id we are using is the same one we report in the otu_dict
                    tax_id = self.data.otu_dict[otu_id].get('^ncbi:taxon') 
                    assert tax_id not in set([0, None]) # every OTU must have a taxon_id for filter blast
                    # we cannot include unmapped taxa in fliter blast.
                if self.downtorank is not None:
                    downtorank_name = None
                    downtorank_id = None
                    if self.config.blast_loc == 'remote':
                        sys.stderr.write("Filtering by taxon not functional for remote ncbi searches yet.")
                        sys.exit(-7)
                        #tax_id = self.ids.get_rank_info_from_web(taxon_name=tax_name)
                        #lineage2ranks = self.ids.otu_rank[tax_id]["rank"]
                        #ncbi = NCBITaxa()
                        #if lineage2ranks == 'unassigned':
                        #    downtorank_id = tax_id
                        #    downtorank_name = tax_name
                        #else:
                        #    for key_rank, val in lineage2ranks.items():
                        #        if val == downtorank:
                        #            downtorank_id = key_rank
                        #            value_d = ncbi.get_taxid_translator([downtorank_id])
                        #            downtorank_name = value_d[int(downtorank_id)]
                    else:
                        downtorank_id = self.ids.ncbi_parser.get_downtorank_id(tax_id, self.downtorank)
                        downtorank_name = self.ids.ncbi_parser.get_name_from_id(downtorank_id)
                    tax_name = downtorank_name
                    tax_id = downtorank_id
                tax_name = tax_name.replace(" ", "_") #TODO we should store the names with underscores in the dicts
                self.ids.spn_to_ncbiid[tax_name] = tax_id
                self.ids.ncbiid_to_spn[tax_id] = tax_name
                tax_id=int(tax_id)
                if tax_id in self.sp_d:
                    self.sp_d[tax_id].append(otu_id)
                else:
                    self.sp_d[tax_id] = [otu_id]
        return self.sp_d

    def make_sp_seq_dict(self):
        """Uses the sp_d to make a dict with species names as key1, key2 is gb_id/sp.name and value is seq

        This is used to select representatives during the filtering step, where it selects how many sequences per
        species to keep in the alignment. It will only contain sp that were not removed in an earlier cycle of the
        program.

        Note: has test, test_sp_seq_d.py

        return: self.sp_seq_d
        """
        debug("make_sp_seq_dict")
        for key in self.sp_d:
            # loop to populate dict. key1 = sp id, key2= gb id, value = seq,
            # number of items in key2 will be filtered according to threshold and already present seq
            seq_d = {}
            for otu_id in self.sp_d[key]:
#                debug("otu id is {}".format(otu_id))
                gb_id = self.data.otu_dict[otu_id].get('^ncbi:accession', otu_id) #this will get the gen bank accession number if it has one, or use the otu_id if not
                try:
                    seq = self.data.aln[otu_id].symbols_as_string()
                except KeyError:
                    assert gb_id in self.new_seqs # if it is not already in the alignment it must be in new seqs
                    seq = self.new_seqs[gb_id]
                seq = seq.replace("-", "")
                seq = seq.replace("?", "")
                seq_d[gb_id] = seq
            self.sp_seq_d[key] = seq_d
        return

    def select_seq_by_local_blast(self, seq_d, fn, count):
        """
        Selects number of sequences from local_blast to fill up sequences to the threshold.
        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program
        how many sequences for the taxon are already available in the aln.

        It will only include species which have a blast score of mean plus/minus sd.
        Uses the information returned by read_local_blast_query()
        to select which sequences shall be added in a filtered run.

        Note: has test,test_select_seq_by_local_blast.py

        :param seq_d: is the value of self.sp_d (= another dict)
        :param fn: refers to a filename to find the local blast file produced before,
                    which needs to be read in by read_local_blast_query() - currently tax_id
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        debug("select_seq_by_local_blast")
        # debug([seq_d, fn, count])
        no_similar_seqs = 0
        seq_blast_score = filter_by_local_blast.read_filter_blast(self.workdir, seq_d, fn)
        random_seq_ofsp = {}

        if seq_blast_score != {}:
            if (self.threshold - count) <= 0:
                debug("already too many samples of sp in aln, skip adding more.")
            # exact amount of seq present which need to be added
            elif len(seq_blast_score.keys()) == (self.threshold - count):
                random_seq_ofsp = seq_blast_score
            elif len(seq_blast_score.keys()) > (
                    self.threshold - count):  # more seq available than need to be added, choose by random
                debug("choose random")
                random_seq_ofsp = random.sample(seq_blast_score.items(), (self.threshold - count))
                random_seq_ofsp = dict(random_seq_ofsp)
            elif len(seq_blast_score.keys()) < (
                    self.threshold - count):  # less seq available than need to be added, just use all
                debug("add all")
                # print(len(seq_blast_score.keys()))
                # print(self.threshold - count)
                random_seq_ofsp = seq_blast_score
        # no similar seq found. think about what to do. was the other seq already present?
        else:
            debug("blast did not find similar seqs")
            if len(seq_d.keys()) > 2 and no_similar_seqs == 0:  # try with different seq to blast
                debug("blast with different seq...")
                # all the next line is from how_many_seq_to_keep()
                blast_seq_id = seq_d.keys()[1]  # seq 1 instead of 0 now
                seq = seq_d[blast_seq_id]
                filter_by_local_blast.write_filterblast_query(self.workdir, blast_seq_id, seq,
                                                              fn=fn)  # blast guy
                blast_db = seq_d.keys()[2:]
                for blast_key in blast_db:
                    seq = seq_d[blast_key]
                    filter_by_local_blast.write_filterblast_db(self.workdir, blast_key, seq,
                                                                  fn=fn)
                # make local blast of sequences
                filter_by_local_blast.run_filter_blast(self.workdir, fn, fn)
                no_similar_seqs = 1
                count_dict = self.count_num_seq(fn)
                seq_present = count_dict["seq_present"]
                if len(seq_d) + seq_present >= self.threshold:
                    self.select_seq_by_local_blast(seq_d, fn, count)
                elif len(seq_d) + seq_present < self.threshold:
                    self.add_all(fn)
            elif len(seq_d.keys()) > 2 and no_similar_seqs == 1:  # also with different seq no result, select random seq!
                if len(seq_d.keys()) == (self.threshold - count):  # exact amount of seq present which need to be added
                    for item in seq_d.keys():
                        random_seq_ofsp[item] = seq_d[item]
                elif len(seq_d.keys()) > (
                        self.threshold - count):  # more seq available than need to be added, choose by random
                    debug("choose random - else")
                    random_seq = random.sample(seq_d.items(), (self.threshold - count))
                    for item in random_seq.keys():
                        random_seq_ofsp[item] = random_seq[item]
                elif len(seq_blast_score.keys()) < (
                        self.threshold - count):  # less seq available than need to be added, just use all
                    debug("add all - else")
                    for item in seq_d.keys():
                        random_seq_ofsp[item] = seq_d[item]
            else:
                seq_id = seq_d.keys()[1]
                seq = seq_d[seq_id]
                random_seq_ofsp[seq_id] = seq
        # debug(random_seq_ofsp)
        if len(random_seq_ofsp) > 0:  # add everything to filtered seq
            for key, val in random_seq_ofsp.items():
                self.filtered_seq[key] = val
        return self.filtered_seq

    def select_seq_by_length(self, taxon_id, count):
        """This is another mode to filter the sequences, if there are more than the threshold.

        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".
        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program how many
        sequences for the taxon are already available in the aln.

        !!! sometimes the only seq in seq_w_maxlen is the original seq,
        then this is the one to be added, but it will be removed,
        later as it is no new seq! thus no new seq for that species is added

        :param taxon_id: key from self.sp_seq_d
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        debug("select_seq_by_length")
        max_len = max(self.sp_seq_d[taxon_id].values())

        seq_w_maxlen = {}
        for key, val in self.sp_seq_d[taxon_id].items():
            for otu_id in self.sp_d[taxon_id]:
                item = self.data.otu_dict[otu_id]
                if '^ncbi:accession' in item and item['^ncbi:accession'] == key:
                    if item['^physcraper:status'].split(' ')[0] != ["added", "deleted", "original", "new"]:
                        if len(val) == len(max_len):
                                seq_w_maxlen[key] = val
        if (self.threshold - count) <= 0:
            debug("already to many samples of sp in aln, skip adding more.")
            random_seq_ofsp = None
        elif len(seq_w_maxlen) == (self.threshold - count):
            random_seq_ofsp = seq_w_maxlen
        elif len(seq_w_maxlen) > (self.threshold - count):
            random_seq_ofsp = random.sample(seq_w_maxlen.items(), (self.threshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        else:
            toselect = range(len(seq_w_maxlen), (self.threshold - count))
            keymax = seq_w_maxlen.keys()
            subdict = {k: v for k, v in self.sp_seq_d[taxon_id].items() if k not in keymax}
            second_len = max(subdict.values())
            seq2len = {}
            for key, val in subdict.items():
                if len(val) == len(second_len):
                    seq2len[key] = val
            random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
            random_seq_ofsp = dict(random_seq_ofsp)
            random_seq_ofsp.update(seq_w_maxlen)
        if (self.threshold - count) >= 1:
            assert random_seq_ofsp is not None
        if random_seq_ofsp is not None:
            for key in random_seq_ofsp.keys():
                self.filtered_seq[key] = random_seq_ofsp[key]

    def add_all(self, key):
        """It adds all seq to filtered_seq dict as the number of sequences present is smaller than the threshold value.

        It is only used, when sequences selection happens via blasting.

        Note: has test, test_add_all.py

        :param key: key of self.sp_d (taxon id)
        :return: self.filtered_seq
        """
        debug('add_all')
        for otu_id in self.sp_d[key]:
            otu_info = self.data.otu_dict[otu_id]
            if '^physcraper:status' in otu_info:
                # debug(otu_id)
                if otu_info['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if otu_info['^physcraper:last_blasted'] == '1800/01/01' \
                            and otu_info['^physcraper:status'] != "original":
                        gb_id = otu_info['^ncbi:accession']
                        seq = self.new_seqs[gb_id]
                        self.filtered_seq[gb_id] = seq
        return self.filtered_seq

    def loop_for_write_blast_files(self, tax_id):
        """This loop is needed to be able to write the local blast files for the filtering step correctly.

        Function returns a filename for the filter blast.

        Note: has test,test_loop_for_blast.py

        :param key: key of self.sp_d (taxon id)
        :return: name of the blast file
        """
        debug("loop_for_write_blast_files")
        aln_otus = set([taxon.label for taxon in self.data.aln])
        query_otu = None
        for otu_id in self.sp_d[tax_id]: 
            otu_info = self.data.otu_dict[otu_id]
            if otu_id in aln_otus:
                query_otu = otu_id #we end up overwriting the query file repeatedly. Might as well just choose one otu and write it once.
            elif '^physcraper:status' in otu_info and otu_info['^physcraper:status'].split(' ')[0] not in self.seq_filter: # these are the new sequences that haven't been filtered out
                assert '^ncbi:accession' in otu_info
                gb_id = otu_info['^ncbi:accession']
                assert len(gb_id.split(".")) == 2 #we should have already gotten rid of any bad ids
                assert gb_id in self.new_seqs.keys()
                seq = self.new_seqs[gb_id]
                filter_by_local_blast.write_filterblast_db(self.workdir, gb_id, seq, fn=tax_id)
            else:
                debug("otu_id {} was not in the alignemnt, but was filtered out due to {}".format(otu_id, otu_info['^physcraper:status']))
        assert query_otu is not None#at least 1 otu must be in the alignment
        debug("for taxon {} will use otu {} for query".format(tax_id, query_otu))
        query_seq = self.data.aln[query_otu]
        filter_by_local_blast.write_filterblast_query(self.workdir, query_otu, query_seq, fn=tax_id)
        return tax_id


    def count_num_seq(self, tax_id):
        """Counts how many sequences there are for a tax_name, excluding sequences that have not been added
        during earlier cycles.

        Function is only used in how_many_sp_to_keep().

        :param tax_id: key from self.sp_seq_d
        :return: dict which contains information of how many seq are already present in aln, how many new seq have been
                found and if the taxon is a new taxon or if seq are already present
        """
        debug("count_num_seq for tax_id {}".format(tax_id))
        seq_added = 0
        original = 0
        new_taxon = True
        query_count = 0
        seq_in_aln = 0
        for otu_id in self.sp_d[tax_id]:
            item = self.data.otu_dict[otu_id]
            aln_otus = set([taxon.label for taxon in self.data.aln])
            if otu_id in aln_otus:
                seq_in_aln += 1
                new_taxon = False
            status = item.get('^physcraper:status')
            assert status is not None
            if status.split(' ')[0] not in self.seq_filter:
#                debug(item['^physcraper:status'])
                item_split = item['^physcraper:status'].split(' ')[0]
                # if item['^physcraper:last_blasted'] != '1800/01/01':
                    # new_taxon = False
                # if item["^physcraper:status"] == "query" or item_split == "new":
                if item["^physcraper:status"] == "query" or item_split == "new" or item_split == "added,":
                    query_count += 1
                if item["^physcraper:status"] == 'added as representative of taxon':
                    seq_added += 1
                    new_taxon = False
                if item_split == "original":
                    original += 1
                    new_taxon = False
        seq_present = seq_added + original
        assert seq_in_aln == seq_present
        # if item_split == "added," or item_split == "original":
        count_dict = {
            "seq_present": seq_added + original,
            "query_count": query_count,
            "new_taxon": new_taxon,
        }
        # if self.config.add_lower_taxa is not True:
        if new_taxon is False:
            assert original != 0 or seq_added != 0, ("count_dict `{}` has more seq added than threshold: 0.".format(count_dict))
        if new_taxon is True:
            assert original == 0, ("count_dict `{}` has more original seq than allowed for new taxon.".format(count_dict))
            assert seq_added == 0, ("count_dict `{}` has more seq added than allowed for new taxon.".format(count_dict))
        # debug([seq_added, original, self.threshold])
        if self.config.add_lower_taxa is not True:
            if original < self.threshold:
                assert seq_added <= self.threshold, ("count_dict `{}` has more seq added than threshold.".format(count_dict))
            elif original > self.threshold:
                sys.stdout.write("already more originals than requested by threshold...\n")
            else:
                assert seq_added + original <= self.threshold, \
                    "seq_added ({}) and original ({}) have more than threshold ({}).".format(seq_added, original, self.threshold)
        return count_dict

    def how_many_sp_to_keep(self, selectby):
        """Uses the sp_seq_d and places the number of sequences according to threshold into the self.filterdseq_dict.

        This is essentially the key function of the Filter-class, it wraps up everything.

        :param selectby: mode of sequence selection, defined in input
        :return: nothing specific, it is the function, that completes the self.filtered_seq, which contains the filtered
                sequences that shall be added.
        """
        debug("how_many_sp_to_keep")
        # self.threshold = threshold
        for tax_id in self.sp_d:
            count_dict = self.count_num_seq(tax_id)
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            new_taxon = count_dict["new_taxon"]
            #debug(count_dict)
            # debug(tax_id)
            if seq_present <= self.threshold:  # add seq to aln
                if seq_present + query_count <= self.threshold:  # to add all stuff to self.filtered_seq[gi_n]
                    debug("taxon {}, has {} seqs present, and we found {}, which totals less than the threshold {}, so we will add all".format(tax_id, seq_present, query_count, self.threshold))
                    self.add_all(tax_id)
                else:  # filter number of sequences
                    debug("taxon {}, has {} seqs present, and we found {}, which totals MORE than the threshold {}, so we will select by {}".format(tax_id, seq_present, query_count, self.threshold,selectby))
                    if tax_id in self.sp_seq_d.keys():
                        if selectby == "length":
                            self.select_seq_by_length(tax_id, seq_present)
                        elif selectby == "blast":
                            if seq_present == 0 and new_taxon is True and query_count >= 1:  # if new taxon
                                # debug("new taxon")
                                blast_seq_id = self.sp_seq_d[tax_id].keys()[0]
                                seq = self.sp_seq_d[tax_id][blast_seq_id]
                                filter_by_local_blast.write_filterblast_query(self.workdir, blast_seq_id, seq,
                                                                              fn=tax_id)  # blast guy
                                blast_db = self.sp_seq_d[tax_id].keys()[1:]
                                for blast_key in blast_db:
                                    seq = self.sp_seq_d[tax_id][blast_key]
                                    filter_by_local_blast.write_filterblast_db(self.workdir, blast_key, seq, fn=tax_id)
                                # make local blast of sequences
                                filter_by_local_blast.run_filter_blast(self.workdir, tax_id, tax_id)
                                if len(self.sp_seq_d[tax_id]) + seq_present >= self.threshold:
                                    self.select_seq_by_local_blast(self.sp_seq_d[tax_id], tax_id, seq_present)
                                elif len(self.sp_seq_d[tax_id]) + seq_present < self.threshold:
                                    self.add_all(tax_id)
                            elif 1 <= seq_present < self.threshold and new_taxon is False and query_count != 0:
                                # species is not new in alignment, make blast with existing seq
                                # debug("old_taxon")
                                # debug([seq_present, query_count])
                                if query_count + seq_present > self.threshold:
                                    taxonfn = self.loop_for_write_blast_files(tax_id)
                                    # debug([tax_id, taxonfn])
                                    if self.downtorank is not None:
                                        taxonfn = tax_id
                                    filter_by_local_blast.run_filter_blast(self.workdir, taxonfn, taxonfn)
                                    self.select_seq_by_local_blast(self.sp_seq_d[tax_id], taxonfn, seq_present)
                                    # debug([tax_id])
                                    # debug(self.filtered_seq)
                                elif query_count + seq_present <= self.threshold:
                                    self.add_all(tax_id)
        # debug(self.filtered_seq)
        return

    def replace_new_seq(self):
        """Function to replace self.new_seqs and self.new_seqs_otu_id with the subset of filtered sequences.

        This is the final step in the FilterBlast class, from here it goes back to PhyScraper.

        :return: subsets of self.new_seqs and self.new_seqs_otu_id
        """
        debug("replace new seq")
        keylist = self.filtered_seq.keys()
        if not self.unpublished:
            keylist = [x for x in keylist if x[:6] != "unpubl"]
        seq_not_added = self.new_seqs.keys()
        reduced_new_seqs_dic = {}
        for gb_id in seq_not_added:
            for key in self.data.otu_dict.keys():
                if '^ncbi:accession' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:accession'] == gb_id:
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        # debug(self.data.otu_dict[key]['^physcraper:status'])
                        if self.data.otu_dict[key]['^physcraper:status'] == "query" \
                                or self.data.otu_dict[key]['^physcraper:status'].split(" ")[0] == 'new':
                            self.data.otu_dict[key]['^physcraper:status'] = 'not added, ' \
                                                                            'there are enough seq per sp in tre'
        for gb_id in keylist:
            added = False
            for key in self.data.otu_dict.keys():
                if "^ncbi:accession" in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]["^ncbi:accession"] == gb_id and added is False:
                        added = True
                        reduced_new_seqs_dic[key] = self.filtered_seq[gb_id]
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added as representative of taxon'

        reduced_new_seqs = {k: self.filtered_seq[k] for k in keylist}
        # debug(reduced_new_seqs_dic)
        with open(self.logfile, "a") as log:
            log.write("{} sequences added after filtering, of {} before filtering\n".format(len(reduced_new_seqs_dic),
                                                                                            len(self.new_seqs_otu_id)))
        self.new_seqs = deepcopy(reduced_new_seqs)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic)
        # set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return




def get_acc_from_blast(query_string):
    """
    Get the accession number from a blast query.
    
    Get acc is more difficult now, as new seqs not always have gi number, then query changes.
    
    :param query_string: string that contains acc and gi from local blast query result
    :return: gb_acc 

    """
    if len(query_string.split("|")) >= 3:
        gb_acc = query_string.split("|")[3]
    else:
        gb_acc = query_string.split("|")[0]
    if len(gb_acc.split(".")) < 2:
        sys.stderr.write("query string {} does not contain a Genbank accession number.".format(query_string))
        exit(-15)
    assert len(gb_acc.split(".")) >= 2, (len(gb_acc.split(".")), gb_acc)
    return gb_acc

def get_gi_from_blast(query_string):
    """
    Get the gi number from a blast query. 
    Get acc is more difficult now, as new seqs not always have gi number, then query changes.

    If not available return None. 

    :param query_string: string that contains acc and gi from local blast query result
    :return: gb_id if available
    """      
    if len(query_string.split("|")) >= 3:
        gb_id = query_string.split("|")[1]
    else:
        return None
    assert len(gb_id.split(".")) < 2, (len(gb_id.split(".")), gb_id)
    assert gb_id.isdigit() is True
    return int(gb_id)



####################
#Funcs below here should be moved to a separate script at some point


def get_tax_info_from_acc(gb_id, data_obj, ids_obj):
    '''takes an accessionumber and returns the ncabi_id and the taxon name'''
#    debug("Getting tax info from acc {}".format(gb_id))
    ncbi_id = None
    tax_name = None
    if gb_id[:6] == "unpubl":  # There may not be ncbi id, because they aren't published
            tax_name = data_obj.gb_dict[gb_id]["^ot:ottTaxonName"]
            ncbi_id = data_obj.gb_dict[gb_id]["^ncbi:taxon"]
            ott_id = data_obj.gb_dict[gb_id]["^ot:ottId"]
            if tax_name is None:
                tax_name = data_obj.gb_dict[gb_id][u'^user:TaxonName']
            if ncbi_id is None: 
                # debug(tax_name.split(" ")[0])
                tax_lin_name = tax_name.split(" ")[0]
                tax_lin_name = tax_lin_name.split("_")[0]
                # debug(tax_lin_name)
                ncbi_id = ids_obj.ncbi_parser.get_id_from_name(tax_lin_name) #TODO What should happen here if the unpublished sequence doesn't have a name that is found?
    elif len(gb_id.split(".")) >= 2:  # used to figure out if gb_id is from Genbank
            if gb_id in data_obj.gb_dict.keys() and "staxids" in data_obj.gb_dict[gb_id].keys():
                tax_name = data_obj.gb_dict[gb_id]["sscinames"]
                ncbi_id = data_obj.gb_dict[gb_id]["staxids"]
            else:  # all web blast results
                tax_name = ids_obj.find_name_acc(acc=gb_id)
                if tax_name is None:
                    sys.stderr.write("no species name returned for {}".format(gb_id))
                ncbi_id = ids_obj.get_ncbiid_from_acc(gb_id)
    else:
        try:
            ncbi_id = ids_obj.get_ncbiid_from_acc(gb_id)
            tax_name = ids_obj.find_name_acc(acc=gb_id)
        except:
            pass
    if ncbi_id == None:
        sys.stderr.write("Failed to get information for sequence with accession number {}".format(gb_id))
    return ncbi_id, tax_name


def get_ncbi_tax_id(handle):
    """Get the taxon ID from ncbi. ONly used for web queries

    :param handle: NCBI read.handle
    :return: ncbi_id
    """
    ncbi_id = None
    gb_list = handle[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
    for item in gb_list:
        if item[u"GBQualifier_name"] == "db_xref":
            if item[u"GBQualifier_value"][:5] == "taxon":
                ncbi_id = int(item[u"GBQualifier_value"][6:])
                break
            else:
                continue
    return ncbi_id


def get_ncbi_tax_name(handle):
    """Get the sp name from ncbi. 
    Could be replaced by direct lookup to ott_ncbi.

    :param handle: NCBI read.handle
    :return: ncbi_spn
    """
    ncbi_sp = None
    gb_list = handle[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
    for item in gb_list:
        if item[u"GBQualifier_name"] == "organism":
            ncbi_sp = str(item[u"GBQualifier_value"])
            ncbi_sp = ncbi_sp.replace(" ", "_")
    return ncbi_sp

