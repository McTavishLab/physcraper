#!/usr/bin/env python
"""Physcraper module"""
import sys
import re
import os
import csv
import subprocess
# import time
import datetime
import glob
import json
# import unicodedata
import configparser
import pickle
# import inspect
import random
# import logging
# import collections
from copy import deepcopy
from ete2 import NCBITaxa
# from urllib2 import URLError
import physcraper.AWSWWW as AWSWWW
# import numpy
from Bio.Blast import NCBIXML  # , NCBIWWW
# from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import Entrez  # , SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from dendropy import Tree, \
    DnaCharacterMatrix, \
    DataSet, \
    datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI, APIWrapper
from peyotl.sugar import tree_of_life, taxomachine  # taxonomy,
from peyotl.nexson_syntax import extract_tree, \
    get_subtree_otus, \
    extract_otu_nexson, \
    PhyloSchema  # extract_tree_nexson, \
# from peyotl.api import APIWrapper

# extension functions
import concat  # is the local concat class
import ncbi_data_parser  # is the ncbi data parser class and associated functions
import local_blast


_DEBUG = 1
_DEBUG_MK = 1
_deep_debug = 0

_VERBOSE = 1


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)
    # with open("debugging.txt", "a") as debugf:
    #     debugf.write("{}\n".format(msg))


def deep_debug(msg):
    """short debugging command
    """
    if _deep_debug == 1:
        print(msg)
    # with open("debugging.txt", "a") as debugf:
    #     debugf.write("{}\n".format(msg))


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


# which python physcraper file do I use?
print("Current --init-- version number: 10-10-2018.0")

debug(os.path.realpath(__file__))


def get_raw_input():
    """Asks for yes or no user input.

    :return: user input
    """
    debug("get raw input")
    is_valid = 0
    while not is_valid:
        try:
            x = raw_input("Please write either 'yes' or 'no': ")
            print(x)
            if x == "yes" or x == "no":
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError, e:
            print ("'%s' is not a valid answer." % e.args[0].split(": ")[1])
    return x


class ConfigObj(object):
    """Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store.

    To build the class the following is needed:
        configfi: a configuration file in a specific format, e.g. to read in self.e_value_thresh.
                    The file needs to have a heading of the format: [blast] and then somewhere below that heading
                    a string e_value_thresh = value

    During the initializing process the following self objects are generated:
        self.e_value_thresh: the defined threshold for the e-value during Blast searches, check out:
                            https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
        self.hitlist_size: the maximum number of sequences retrieved by a single blast search
        self.seq_len_perc: value from 0 to 1. Defines how much shorter new seq can be compared to input
        # self.get_ncbi_taxonomy: Path to sh file doing something...
        # self.ncbi_dmp: path to file that has gi numbers and the corresponding ncbi tax id's
        self.phylesystem_loc: Default is to run on remote, github phylesystem, can be set to 'local'
                              to access files from local clone
        self.ott_ncbi: path to file containing OTT id, ncbi and taxon name
        self.id_pickle: path to pickle file 
        self.email: email address used for blast queries
        self.blast_loc: defines which blasting method to use, either web-query (=remote) or from a local
                        blast database (=local)
        self.num_threads number of cores to be used during a run
        self.gifilename: True or False, defines if blast results shall be shared across runs
        self.url_base: if blastloc == remote, it defines the url for the blast queries.
                        if blastloc == local: url_base = None

        optional self.objects:
            self.blastdb: if blastloc == local, this defines the path to the local blast database

        self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
        self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
        self.unmapped: remove/keep - used for OToL original tips that can not be assigned to a taxon.
                        remove - tips will be removed/ keep - tip will be assigned to mrca id
    """

    def __init__(self, configfi):
        if _DEBUG:
            sys.stdout.write("Building config object\n")
        debug(configfi)
        debug(os.path.isfile(configfi))
        assert os.path.isfile(configfi)
        config = configparser.ConfigParser()
        config.read(configfi)
        self.e_value_thresh = config['blast']['e_value_thresh']
        assert is_number(self.e_value_thresh)
        self.hitlist_size = int(config['blast']['hitlist_size'])
        self.seq_len_perc = float(config['physcraper']['seq_len_perc'])
        assert 0 < self.seq_len_perc < 1
        # self.get_ncbi_taxonomy = config['taxonomy']['get_ncbi_taxonomy']  # TODODELETE: Was dropped in mov from subprocess I think/...
        # assert os.path.isfile(self.get_ncbi_taxonomy)
        # self.ncbi_dmp = config['taxonomy']['ncbi_dmp']  # TODODELETE: Currently not used. Was used to get tax_id from gi via subprocess. 10GB!
        # gi_id to taxid (according to GenBank it's not updated since 2016, even though the files seems to be newer)
        # if not os.path.isfile(self.ncbi_dmp):
        #     os.system("rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz {}.gz".format(self.ncbi_dmp))
        #     os.system("gunzip taxonomy/gi_taxid_nucl.dmp.gz")
        #     self.ncbi_dmp = "taxonomy/gi_taxid_nucl.dmp.gz"
        self.phylesystem_loc = config['phylesystem']['location']
        assert (self.phylesystem_loc in ['local', 'api'])  # default is api, but can run on local version of OpenTree datastore
        self.ott_ncbi = config['taxonomy']['ott_ncbi']
        assert os.path.isfile(self.ott_ncbi)
        self.id_pickle = os.path.abspath(config['taxonomy']['id_pickle'])  # rewrites the relative path as an absolute path so that it behaves when changing dirs
        self.email = config['blast']['Entrez.email']
        assert '@' in self.email
        self.blast_loc = config['blast']['location']
        self.num_threads = config['blast'].get('num_threads')
        assert self.blast_loc in ['local', 'remote']
        if self.blast_loc == 'local':
            self.blastdb = config['blast']['localblastdb']
            self.url_base = None
        if self.blast_loc == 'remote':
            self.url_base = config['blast'].get('url_base')
        self.gifilename = config['blast'].get('gifilename', False)
        if self.gifilename is not False:
            if self.gifilename == "True" or self.gifilename == "true":
                self.gifilename = True
            else:
                self.gifilename = False
        self.ncbi_parser_nodes_fn = config['ncbi_parser']["nodes_fn"]
        self.ncbi_parser_names_fn = config['ncbi_parser']["names_fn"]
        # print("check file status")
        self._download_ncbi_parser()
        self._download_localblastdb()
        self.unmapped = config['blast']['unmapped']
        assert self.unmapped in ['remove', 'keep']
        if _DEBUG:
            sys.stdout.write("{}\n".format(self.email))
            if self.blast_loc == 'remote':
                sys.stdout.write("url base = {}\n".format(self.url_base))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == 'local':
                sys.stdout.write("local blast db {}\n".format(self.blastdb))


    def _download_localblastdb(self):
        """Check if files are present and if they are uptodate.
        If not files will be downloaded. 
        """
        if self.blast_loc == 'local':
            if not os.path.isfile("{}/nt.01.nhr".format(self.blastdb)):
                print("Do you want to download the blast nt databases from ncbi? Note: "
                      "This is a US government website! You agree to their terms")
                x = get_raw_input()
                if x == "yes":
                    os.system("rsync -av ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*" 
                              "{}/".format(self.blastdb))
                    os.system("rsync -av ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz" 
                              "{}/".format(self.blastdb))
                    cwd = os.getcwd()
                    os.chdir(self.blastdb)
                    os.system("update_blastdb nt")
                    os.system("cat *.tar.gz | tar -xvzf - -i")
                    os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                    os.chdir(cwd)
                elif x == "no":
                    print("You did not agree to download data from ncbi. Programm will default to blast web-queries.")
                    print("This is slow and crashes regularly!")
                    self.blast_loc = 'remote'
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime("{}/nt.01.nhr".format(self.blastdb))
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days    
                # debug([download_date, today, time_passed])
                if time_passed >= 60: 
                    print("Your databases might not be uptodate anymore. You downloaded them {} days ago. "
                          "Do you want to update the blast databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms".format(time_passed))
                    x = get_raw_input()
                    if x == "yes":
                        os.system('update_blastdb nt') 
                        os.system("perl update_blastdb.pl taxdb")
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type 'yes' or 'no'!")

    def _download_ncbi_parser(self):
        """Check if files are present and if they are uptodate.
        If not files will be downloaded. 
        """
        if self.blast_loc == 'local':
            if not os.path.isfile(self.ncbi_parser_nodes_fn):
                print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms")
                x = get_raw_input()
                if x == "yes":
                    os.system("rsync -av ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" 
                              "./tests/data/taxdump.tar.gz")
                    os.system("gunzip - cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                elif x == "no":
                    print("You did not agree to download data from ncbi. Programm will default to blast web-queries.")
                    print("This is slow and crashes regularly!")
                    self.blast_loc = 'remote'
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime(self.ncbi_parser_nodes_fn)
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days    
                # debug([download_date, today, time_passed])
                if time_passed >= 60: 
                    print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms")
                    x = get_raw_input()
                    if x == "yes":
                        os.system("rsync -av ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" 
                                  "./tests/data/taxdump.tar.gz")
                        os.system("gunzip - cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type yes or no!")


# ATT is a dumb acronym for Alignment Tree Taxa object
def get_dataset_from_treebase(study_id,
                              phylesystem_loc='api'):
    """Function is used to get the aln from treebase, for a tree that OpenTree has the mapped tree.
    """
    try:
        nexson = get_nexson(study_id, phylesystem_loc)
    except:  # TODO: seems to be an http error. Did not fgure out how to handle them (requests.exceptions.HTTPError)
        sys.stderr.write("couldn't find study id {} in phylesystem location {}\n".format(study_id, phylesystem_loc))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treebase record associated with study ")
        sys.exit()
    else:
        tb_id = treebase_url.split(':S')[1]
        url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexus".format(tb_id)
        # url = "http://treebase.org/treebase-web/phylows/study/TB2:S{}?format=nexml".format(tb_id)
        if _DEBUG:
            sys.stderr.write(url + "\n")
        dna = DataSet.get(url=url,
                          schema="nexml")
        return dna


def generate_ATT_from_phylesystem(aln,
                                  workdir,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api',
                                  ingroup_mrca=None):
    """gathers together tree, alignment, and study info - forces names to otu_ids.
    Outputs AlignTreeTax object.

    Input can be either a study ID and tree ID from OpenTree  
    # TODO: Is there another way to get those data, than actually going to the OToL website, 
    or can I query the internet for a desired part of the tree and corresponding studies?TODO MK: write script to easily get those information: Yes it can, use peyotl find_trees
     
     According to code it cannot be either, but must be both

    Alignemnt need to be a Dendropy DNA character matrix!

    :param aln: dendropy alignment object
    :param workdir: path to working directory
    :param study_id: OToL study id of the corresponding phylogeny which shall be updated
    :param tree_id: OToL corresponding tree ID as some studies have several phylogenies
    :param phylesystem_loc: access the github version of the OpenTree data store, or a local clone
    :param ingroup_mrca: OToL identifier of the mrca of the clade that shall be updated (can be subset of the phylogeny)
    :return: object of class ATT
    """
    # TODO CHECK ARGS
    assert isinstance(aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore UGH
    nexson = get_nexson(study_id, phylesystem_loc)
    ott_ids = get_subtree_otus(nexson,
                               tree_id=tree_id,
                               subtree_id="ingroup",
                               return_format="ottid")
    if ingroup_mrca:
        if type(ingroup_mrca) == list:
            ott_ids = set(ingroup_mrca)
            # debug(ott_ids)

            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(ingroup_mrca)
    else:
        ott_mrca = get_mrca_ott(ott_ids)
    newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
    newick = newick.replace(" ", "_")  # UGH Very heavy handed, need to make sure happens on alignement side as well.
    tre = Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    # print("get_subtree_otus")
    otus = get_subtree_otus(nexson, tree_id=tree_id) #this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    # print(otus)
    otu_dict = {}
    orig_lab_to_otu = {}
    treed_taxa = {}
    for otu_id in otus:
        # print(otu_id)
        # print(extract_otu_nexson(nexson, otu_id))
        otu_dict[otu_id] = extract_otu_nexson(nexson, otu_id)[otu_id]
        otu_dict[otu_id]['^physcraper:status'] = "original"
        otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"
        orig = otu_dict[otu_id].get(u'^ot:originalLabel').replace(" ", "_")
        orig_lab_to_otu[orig] = otu_id
        treed_taxa[orig] = otu_dict[otu_id].get(u'^ot:ottId')
    for tax in aln.taxon_namespace:
        try:
            tax.label = orig_lab_to_otu[tax.label].encode('ascii')
        except KeyError:
            sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignment. "
                             "This may indicate a mismatch between tree and alignment\n".format(tax.label))
    # need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    workdir = os.path.abspath(workdir)
    # print(treed_taxa)
    # print(some)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)
    # newick should be bare, but alignment should be DNACharacterMatrix


# def convert(data):
#     """convert json 2.7 as string problem"""
#     # TODO MK: seems not to be used anymore. Was used to convert json strings, as they had a weird structure.
#     if isinstance(data, basestring):
#         return str(data)
#     elif isinstance(data, collections.Mapping):
#         return dict(map(convert, data.iteritems()))
#     elif isinstance(data, collections.Iterable):
#         return type(data)(map(convert, data))
#     else:
#         return data


def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            treefile,
                            otu_json,
                            schema_trf,
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem.

    If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.
    otu_json should encode the taxon names for each tip.

    Note: has test -> owndata.py

    :param seqaln: path to sequence alignment
    :param mattype: string containing format of sequence alignment
    :param workdir: path to working directory
    :param treefile: path to phylogeny
    :param otu_json: path to jsonfile containing the translation of tipnames to taxon names
    :param schema_trf: string defining the format of the input phylogeny
    :param ingroup_mrca: optional - OToL ID of the mrca of the clade of interest
    :return: object of class ATT
    """

    # replace ? in seqaln with - : papara handles them as different characters
    with open(seqaln, 'r') as fin:
        filedata = fin.read()
    filedata = filedata.replace('?', '-')
    # Write the file out again
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    new_seq_file = "{}/replaced_inputaln.fasta".format(workdir)
    with open("{}/replaced_inputaln.fasta".format(workdir), 'w') as aln_file:
        aln_file.write(filedata)
    # use replaced aln as input
    aln = DnaCharacterMatrix.get(path=new_seq_file, schema=mattype)
    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore UGH
    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    assert tre.taxon_namespace is aln.taxon_namespace
    otu_newick = tre.as_string(schema=schema_trf)
    otu_dict = json.load(open(otu_json, "r"))
    debug("get mrca")
    debug(ingroup_mrca)
    # debug(some)
    if ingroup_mrca:
        if type(ingroup_mrca) == list:
            ott_ids = set(ingroup_mrca)
            # debug(ott_ids)

            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(ingroup_mrca)
    else:
        ott_ids = [otu_dict[otu].get(u'^ot:ottId', ) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)


def standardize_label(item):
    """Make sure that the tipnames are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

    Note: has test -> test_edit_dict_key.py

    :param item: original tipname
    :return: tipname in unicode
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
    """
    try:
        res = taxomachine.TNRS(spp_name)['results'][0]
    except IndexError:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0
    if res['matches'][0]['is_approximate_match'] == 1:
        sys.stderr.write("""exact match to taxon {} not found in open tree taxonomy.
                          Check spelling. Maybe {}?""".format(spp_name, res['matches'][0][u'ot:ottTaxonName']))
        return 0
    if res['matches'][0]['is_approximate_match'] == 0:
        ottid = res['matches'][0]['taxon'][u'ott_id']
        ottname = res['matches'][0]['taxon'][u'unique_name']
        ncbi_id = None
        for source in res['matches'][0]['taxon'][u'tax_sources']:
            if source.startswith('ncbi'):
                ncbi_id = source.split(":")[1]
        return ottid, ottname, ncbi_id
    else:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
        return 0


def OtuJsonDict(id_to_spn, id_dict):
    """Make otu json dict, which is also produced within the openTreeLife-query
    reads input file into the var sp_info_dict, translates using an IdDict object
    using web to call Open tree, then ncbi if not found.

    This function is used, if files that shall be updated are not part of the OpenTreeofLife project.
    It reads in the file that contains the tipnames and the corresponding species names.
    It then tries to get the different identifier from the OToL project or if not from ncbi.

    :param id_to_spn: user file, that contains tipname and corresponding sp name for input files.
    :param id_dict: uses the id_dict generates earlier
    :return: dictionary with key: "otu_tiplabel" and value is another dict with the keys '^ncbi:taxon',
                                                    '^ot:ottTaxonName', '^ot:ottId', '^ot:originalLabel',
                                                    '^user:TaxonName', '^physcraper:status', '^physcraper:last_blasted'
    """
    sys.stdout.write("Set up OtuJsonDict \n")
    sp_info_dict = {}
    nosp = []
    with open(id_to_spn, mode='r') as infile:
        for lin in infile:
            ottid, ottname, ncbiid = None, None, None
            tipname, species = lin.strip().split(',')
            # debug(tipname)
            clean_lab = standardize_label(tipname)
            assert clean_lab not in sp_info_dict
            otu_id = "otu{}".format(clean_lab)
            spn = species.replace("_", " ")
            debug(spn)

            info = get_ott_taxon_info(spn)
            if info:
                ottid, ottname, ncbiid = info
            if not info:
                # if _DEBUG:
                #     sys.stdout.write("match to taxon {} not found in open tree taxonomy. "
                #                      "Trying NCBI next.\n".format(spn))
                ncbi = NCBITaxa()
                name2taxid = ncbi.get_name_translator([spn])
                # debug(name2taxid)
                if len(name2taxid.items()) >= 1:
                    # if _DEBUG:
                    #     sys.stdout.write("found taxon {} in ncbi".format(spn))
                    ncbiid = name2taxid.items()[0][1][0]
                else:
                    sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                     "Proceeding without taxon info\n".format(spn))
                    nosp.append(spn)
            sp_info_dict[otu_id] = {'^ncbi:taxon': ncbiid, '^ot:ottTaxonName': ottname, '^ot:ottId': ottid,
                                    '^ot:originalLabel': tipname, '^user:TaxonName': species,
                                    '^physcraper:status': 'original', '^physcraper:last_blasted': "1900/01/01"}
    print(nosp)
    # print(some)
    return sp_info_dict


class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
    Hypothetically, all the keys in the  otu_dict should be clean.

    To build the class the following is needed:
        newick: dendropy.tre.as_string(schema=schema_trf) object
        otu_dict: json file including the otu_dict information generated earlier
        alignment: dendropy DNACharacterMatrix object
        ingroup_mrca: OToL identifier of the group of interest, either subclade as defined by user or of
                        all tiplabels in the phylogeny
        workdir: the path to the corresponding working directory
        schema: optional argument to define tre file schema, if different from "newick"

    During the initializing process the following self objects are generated:
        self.aln: contains the alignment and which will be updated during the run
        self.tre: contains the phylogeny, which will be updated during the run
        self.otu_dict: dictionary with taxon information and physcraper relevant stuff
                key: a unique identifier (otu plus either "tiplabel of phylogeny" or for newly found sequences
                    PS_number.
                value: dictionary:
                        keys: values
                        '^ncbi:gi': GenBank identifier - deprecated by Genbank - only older sequences will have it
                        '^ncbi:accession': Genbanks accession number
                        '^ncbi:title': title of Genbank sequence submission
                        '^ncbi:taxon': ncbi taxon identifier
                        '^ot:ottId': OToL taxon identifier
                        '^physcraper:status': contains information if it was 'original', 'queried', 'removed',
                                            'added during filtering process'
                        '^ot:ottTaxonName': OToL taxon name
                        '^physcraper:last_blasted': contains the date when the sequence was blasted.
                                                    If the year is different from the 20th century, it tells us
                                                    something about the initial status:
                                                     - 1800 = never blasted, not yet considered to be added
                                                     - 1900 = never blasted and not added - see status for more info
                                                     - this century = blasted and added.
                        '^user:TaxonName': optional, user given label from OtuJsonDict
                        "^ot:originalLabel" optional, user given tip label of phylogeny
        self.ps_otu: iterator for new otu IDs, is used as key for self.otu_dict
        self.workdir: contains the path to the working directory, if folder does not exists it is generated.
        self.ott_mrca: OToL taxon Id for the most recent common ancestor of the ingroup
        self.orig_seqlen: list of the original sequence length of the input data
        self.gi_dict: dictionary, that has all information from sequences found during the blasting.
                    key: GenBank sequence identifier
                    value: dictionary, content depends on blast option, differs between webquery and local blast queries
                        key - value pairs for local blast:
                            '^ncbi:gi': GenBank sequence identifier
                            'accession': GenBank accession number
                            'staxids': Taxon identifier
                            'sscinames': Taxon species name
                            'pident': Blast  percentage of identical matches
                            'evalue': Blast e-value
                            'bitscore': Blast bitscore, used for FilterBlast
                            'sseq': corresponding sequence
                            'title': title of Genbank sequence submission
                        key - value for web-query:
                            'accession':Genbank accession number
                            'length': length of sequence
                            'title': string combination of hit_id and hit_def
                            'hit_id': string combination of gi id and accession number
                            'hsps': Bio.Blast.Record.HSP object
                            'hit_def': title from GenBank sequence
                        optional key - value pairs for unpublished option:
                            'localID': local sequence identifier
        # self.orig_aln: original input alignment
        # self.orig_newick: original input phylogeny
        self._reconciled = True/False,
        self.unpubl_otu_json: optional, will contain the OTU-dict for unpublished data, if that option is used

    Following functions are called during the init-process:
        self._reconcile_names(): ## TODO: removes taxa, that are not found in both, the phylogeny and the aln and changes their names????

    The physcraper class is then updating: self.aln, self.tre and self.otu_dict, self.ps_otu, self.gi_dict
    """

    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None, taxon_namespace=None):
        # TODO add assertions that inputs are correct type!!!
        print("build ATT class")
        self.aln = alignment
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
        # debug(self.tre.taxon_namespace)
        assert (self.tre.taxon_namespace is self.aln.taxon_namespace)
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
        # self.tre = Tree.get(data=newick,
        #                     schema="newick",
        #                     preserve_underscores=True,
        #                     taxon_namespace=self.aln.taxon_namespace)
        assert isinstance(otu_dict, dict)
        self.otu_dict = otu_dict
        self.ps_otu = 1  # iterator for new otu IDs
        self._reconcile_names()
        self.workdir = os.path.abspath(workdir)  # TODO: is this where the workdir should live?
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca)
        self.ott_mrca = ingroup_mrca  # TODO: we only use .ott_mrca to infer mrca_ncbi. Why not using the ncbi one directly?
        self.orig_seqlen = []  # FIXME
        self.gi_dict = {}  # has all info about new blast seq  TODODELTE (should maybe go anyhow due to gi switch?): Should this not be part of physcraper class instead? it has all blast information. Blast is not part of this class.
        # self.orig_aln = alignment  # TODODELETE: we never do anything with it.
        # self.orig_newick = newick  # TODODELETE: we never do anything with it.
        self._reconciled = False  # TODO: for what do we want to use it? .... it was checking to see if name reconcilation has ahappened yet. Should get flipped to true when done. MK: Yes, but we never do anything with the information
        self.unpubl_otu_json = None
    
    def _reconcile_names(self):
        """This checks that the tree "original labels" from phylsystem
        align with those found in the alignment. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!
        Taxa that are only found in the tree, or only in the alignment are deleted.
        """
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
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_tax)
            if taxon in aln_tax:
                self.aln.remove_sequences(prune)
            if taxon in treed_tax:
                self.tre.prune_taxa(prune)
        for tax in prune:
            #potentially slow at large number of taxa and large numbers to be pruned
            found = 0
            for otu in self.otu_dict:
                if self.otu_dict[otu][u'^ot:originalLabel'] == tax.label:
                    self.otu_dict[otu]['^physcraper:status'] = "deleted in name reconciliation"
                    found = 1
            if found == 0:
                sys.stderr.write("lost taxon {} in name reconcilliation".format(tax.label))
            self.aln.taxon_namespace.remove_taxon(tax)
        assert (self.aln.taxon_namespace == self.tre.taxon_namespace)
        for tax in self.aln.taxon_namespace:
            # debug(self.otu_dict.keys())
            # debug(tax.label)
            if tax.label in self.otu_dict.keys():
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ''
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for otu in self.otu_dict:
                    # debug(self.otu_dict[otu].get('^ot:originalLabel'))
                    if self.otu_dict[otu].get('^ot:originalLabel') == tax.label or self.otu_dict[otu].get('^ot:originalLabel') == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))
                # debug(some)
                # assert tax.label in self.otu_dict
        # debug(some)
    # TODO - make sure all taxon labels are unique OTU ids.

    def prune_short(self, min_seqlen_perc=0.75):
        """Prunes sequences from alignment if they are shorter than 75%, or if tip is only present in tre.

        Sometimes in the de-concatenating of the original alignment
        taxa with no sequence are generated or in general if certain sequences are really short. In other cases there might be too many tips in the tre
        This gets rid of those from both the tre and the alignment.

        has test: test_prune_short.py

        :param min_seqlen_perc: minimum length of seq
        :return: prunes aln and tre
        """
        # TODO: is this not half the part of _reconcile_names? 
        print(sum(self.orig_seqlen), len(self.orig_seqlen))
        if sum(self.orig_seqlen) != 0:
            avg_seqlen = sum(self.orig_seqlen) / len(self.orig_seqlen)
            seq_len_cutoff = avg_seqlen * min_seqlen_perc
        else:
            for tax, seq in self.aln.items():
                seqlen = len(self.aln[tax].symbols_as_string())
                break
            seq_len_cutoff = seqlen * min_seqlen_perc
        prune = []
        aln_ids = set()
        for tax, seq in self.aln.items():
            aln_ids.add(tax.label)

            if len(seq.symbols_as_string().translate(None, "-?")) <= seq_len_cutoff:
                prune.append(tax)
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
        print(prune)
        if prune:
            # debug(prune)

            # self.tre.prune_taxa(prune)
            # self.aln.taxon_namespace.remove_taxon_label(tax)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short "
                     "step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                # self.aln.remove_sequences(prune)
                # self.tre.prune_taxa_with_labels(prune)  # sometimes it does not delete it with the statement before. Tried to figure out why, have no clue yet.
                self.remove_taxa_aln_tre(tax.label)
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        debug(self.aln.taxon_namespace)
        for tax in prune:
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in prune short"
            # self.aln.taxon_namespace.remove_taxon(tax.label)  # remove_taxon:raises no error, remove_taxon_label: raises error
            # self.tre.taxon_namespace.remove_taxon(tax.label)  # raises error if not found, instead of remove_taxon

            # self.remove_taxa_aln_tre(tax.label)
        debug([self.aln.taxon_namespace, len(self.aln.taxon_namespace)])
        debug([self.tre.taxon_namespace, len(self.tre.taxon_namespace)])
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        # debug([treed_taxa, len(treed_taxa)])
        # debug([aln_ids, len(aln_ids)])
        debug([item for item in treed_taxa if item not in aln_ids])
        debug([item for item in aln_ids if item not in treed_taxa])
        assert treed_taxa.issubset(aln_ids)

        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        # self.reconcile()
        # for key in  self.otu_dict.keys():
        #      if key not in aln_ids:
        #           sys.stderr.write("{} was in otu dict but not alignment. it should be in new seqs...\n".format(key)
        self.trim()
        self._reconciled = 1

    # def reconcile(self, seq_len_perc=0.75):
    #     #TODO RENAME This is pruning short seqs - not reconciling!
    #     """Removes sequences from data if seq is to short.
    #
    #     all missing data seqs are sneaking in, but from where?!
    #
    #     not only used in the beginning...is used to remove sequences that are shorter than 75%"""
    #     # TODO: How is this different from prune_short/ reconcile names?
    #     if prune:
    #         # debug(prune)
    #     for tax in prune:
    #         # debug(tax)
    #         # debug(tax.label)
    #         # debug(self.otu_dict[tax.label])
    #         self.otu_dict[tax.label]['^physcraper:status'] = "deleted in reconcile"
    #         # TODO: line above is unnecessary? get's overwritten in next line by remove_taxa_aln_tre
    #         self.remove_taxa_aln_tre(tax.label)
    #     aln_ids = set()
    #     for tax in self.aln:
    #         aln_ids.add(tax.label)
    #     # debug(len(self.otu_dict.keys()))
    #     debug(len(aln_ids))
    #     debug([item for item in self.otu_dict.keys() if item not in aln_ids])
    #     debug([item for item in aln_ids if item not in self.otu_dict.keys()])
    #     assert aln_ids.issubset(self.otu_dict.keys())
    #     treed_taxa = set()
    #     orphaned_leafs = set()
    #     # here leaf_nodes have taxa that were dropped before. Why do we have this anyways?
    #     for leaf in self.tre.leaf_nodes():
    #         treed_taxa.add(leaf.taxon.label)
    #         if leaf.taxon.label not in aln_ids:
    #             debug(self.otu_dict.keys())
    #             debug(leaf.taxon.label)
    #             self.otu_dict[leaf.taxon.label]['^physcraper:status'] = "deleted due to presence in tree but not aln ?!"
    #             orphaned_leafs.add(leaf)
    #             # TODO figure out why sometimes one of them works and not the other and vice versa
    #             self.tre.prune_taxa([leaf])
    #             # self.tre.prune_taxa_with_labels([leaf.taxon.label])
    #             self.tre.prune_taxa_with_labels([leaf.taxon.label])
    #             self.tre.prune_taxa_with_labels([leaf])
    #             treed_taxa.remove(leaf.taxon.label)
    #             # debug(self.otu_taxonlabel_problem.keys())
    #         # else:
    #         #     treed_taxa.add(leaf.taxon.label)
    #     # debug('treed_taxa')
    #     # debug(treed_taxa)
    #     # debug('aln_ids')
    #     # debug(aln_ids)
    #     # debug([item for item in treed_taxa if item not in aln_ids])
    #     # #debug([item for item in aln_ids if item not in treed_taxa])
    #     # debug(self.tre.taxon_namespace) # otu is gone from namespace, but in treed
    #     # debug(self.tre.as_string(schema="newick")) # otu is in tre, thus  not removed in remove_taxa_aln_tre
    #     assert treed_taxa.issubset(aln_ids)
    #     # for key in  self.otu_dict.keys():
    #     #      if key not in aln_ids:
    #     #           sys.stderr.write("{} was in otu dict but not alignment. it should be in new seqs...\n".format(key)
    #     self.trim()
    #     self._reconciled = 1

    def trim(self, taxon_missingness=0.75):
        """ It removes bases at the start and end of alignments, if they are represented by less than 75%
        of the sequences in the alignment.
        Ensures, that not whole chromosomes get dragged in. It's cutting the ends of long sequences.

        Used in prune_short()
        has test: test_trim.py

        :taxon_missingness: defines how much longer sequence can be in percent
        """
        debug('in trim')
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
            counts = {'?': 0, '-': 0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ['?', '-']:
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
        aln_ids = set()
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        # for tax in self.aln:
            aln_ids.add(taxon.label)
        assert aln_ids.issubset(self.otu_dict.keys())
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        # debug(treed_taxa)
        for leaf in self.tre.leaf_nodes():
            if leaf.taxon not in aln_ids:
                # debug("leaf.taxon not present in aln_ids")
                # debug(leaf.taxon)
                self.tre.prune_taxa([leaf])
                self.tre.prune_taxa_with_labels([leaf.taxon])
                self.tre.prune_taxa_with_labels([leaf])
                # scrape.data.tre.taxon_namespace.remove_taxon_label(leaf.taxon.label)
                treed_taxa.remove(leaf.taxon)
        assert treed_taxa.issubset(aln_ids)
        if _VERBOSE:
            sys.stdout.write("trimmed alignment ends to < {} missing taxa, "
                             "start {}, stop {}\n".format(taxon_missingness, start, stop))
        return

    def check_double_gi(self):
        """for deeper debugging to make sure nothing is added twice"""
        debug("check_double_gi")
        ncbigi_list = []
        for key, val in self.otu_dict.items():
            if '^ncbi:gi' in val:
                gi_otu_dict = val["^ncbi:gi"]
                ncbigi_list.append(str(gi_otu_dict))
        return ncbigi_list

    def add_otu(self, gi_id, ids_obj):
        """ Generates an otu_id for new sequences and adds them into self.otu_dict.
        Needs to be passed an IdDict to do the mapping

        :param gi_id: the Genbank identifier
        :param ids_obj: needs to IDs class to have access to the taxonomic information
        :return: the unique otu_id - the key from self.otu_dict of the corresponding sequence
        """
        # debug("add_otu function")
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        # debug(gi_id)
        ncbi_id = None
        tax_name = None
        ott_id = None
        if type(gi_id) == int:
            # debug("gi_id is int")
            if gi_id in self.gi_dict.keys() and 'staxids' in self.gi_dict[gi_id].keys():
                tax_name = self.gi_dict[gi_id]['sscinames']
                ncbi_id = self.gi_dict[gi_id]['staxids']
                # debug(ncbi_id)
            else:
                tax_name = ids_obj.find_name(gi=gi_id)
                if tax_name is None:
                    sys.stderr.write("no species name returned for {}".format(gi_id))
                ncbi_id = ids_obj.map_gi_ncbi(gi_id)
        elif gi_id[:6] == "unpubl":
            # debug(self.gi_dict.keys())
            # debug(self.gi_dict[gi_id].keys())
            tax_name = self.gi_dict[gi_id]['^ot:ottTaxonName']
            ncbi_id = self.gi_dict[gi_id]['^ncbi:taxon']
            ott_id = self.gi_dict[gi_id]['^ot:ottId']
        else:
            sys.stderr.write("Something is wrong, I cannot add a new seq which has no gi_id or is not unpublished.")

        if _deep_debug == 1:
            currentgilist = self.check_double_gi()
            if gi_id in currentgilist:
                exit(-1)

        if ncbi_id is None:
            debug("ncbi_id is none")

            if ids_obj.otu_rank is not None:
                ncbi_id = ids_obj.otu_rank[tax_name]["taxon id"]
            else:
                ncbi_id = ids_obj.ncbi_parser.get_id_from_name(tax_name)
            # if type(gi_id) == int:
            #     print("add id to self")
            ids_obj.gi_ncbi_dict[gi_id] = ncbi_id
            ids_obj.ncbiid_to_spn[ncbi_id] = tax_name
            ids_obj.spn_to_ncbiid[tax_name] = ncbi_id

        if ncbi_id in ids_obj.ncbi_to_ott.keys():
            # ncbi_id = int(ids_obj.map_gi_ncbi(gi_id))
            ott_id = int(ids_obj.ncbi_to_ott[ncbi_id])
        if ott_id is None:
            ott_id = "OTT_{}".format(self.ps_otu)
            self.ps_otu += 1
            # tax_name = str(ids_obj.ott_to_name[ott]).replace(" ", "_")  # seems to be unused
        if otu_id in self.otu_dict.keys():
            ott_name = ids_obj.ott_to_name.get(ott_id)
        else:
            ott_name = tax_name
        # debug("otu_id")
        # debug(otu_id)
        # debug("ncbi_id")
        # debug(ncbi_id)
        # debug(self.gi_dict[gi_id].keys())

        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi_id
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi_id]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi_id]['title']
        self.otu_dict[otu_id]['^ncbi:taxon'] = ncbi_id
        self.otu_dict[otu_id]['^ot:ottId'] = ott_id
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = ott_name
        # last_blasted date infos: 1800 = never blasted; 1900 = blasted 1x, not added; this century = blasted and added
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"
        if type(gi_id) != int:
            # debug(gi_id)
            # self.otu_dict[otu_id]['^user:TaxonName']
            self.otu_dict[otu_id]['^physcraper:status'] = "local seq"
            self.otu_dict[otu_id]["^ot:originalLabel"] = self.gi_dict[gi_id]['localID']
        if _DEBUG >= 2:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi_id, otu_id))
        # debug(otu_id)
        return otu_id

    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """This writes out needed files for papara (except query sequences).
        Papara is finicky about trees and needs phylip format for the alignment.

        Is only used within func align_query_seqs."""
        # TODO: CAN I even evaluate things in the function definitions?
        # TODO: names for tree and aln files should not be changed, as they are hardcoded in align_query_seqs().
        debug('write papara files')
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
        """Outputs both the streaming files and a ditechecked"""
        # TODO: ditchecked?
        # First write rich annotation json file with everything needed for later?
        debug("write_files")
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema, unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)

    def write_labelled(self, label, treepath=None, alnpath=None, norepeats=True, gi_id=False):
        """output tree and alignment with human readable labels
        Jumps through a bunch of hoops to make labels unique.

        NOT MEMORY EFFICIENT AT ALL

        Has different options available for different desired outputs

        :param label: which information shall be displayed in labelled files: possible options:
                    '^ot:ottTaxonName', '^user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"
        :param treepath: optional: full filenname (including path) for phylogeny
        :param alnpath:  optional: full filenname (including path) for alignment
        :param norepeats: optional: if there shall be no duplicate names in the labelled output files
        :param gi_id: optional, to supplement tiplabel with corresponding GenBank sequence identifier
        :return: writes out labelled phylogeny and alignment to file
        """
        debug("write labelled files")
        if treepath is None:
            treepath = "{}/{}".format(self.workdir, 'labelled.tre')
        if alnpath is None:
            alnpath = "{}/{}".format(self.workdir, 'labelled.aln')
        assert label in ['^ot:ottTaxonName', '^user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
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
            # debug(taxon)
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label is None:
                if self.otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_{}".format(self.otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"),
                                                            self.otu_dict[taxon.label].get('^ot:ottTaxonName', "unk"))
            new_label = str(new_label).replace(' ', '_')
            if gi_id:
                gi_id = self.otu_dict[taxon.label].get('^ncbi:gi')
                if gi_id is None:
                    gi_id = self.otu_dict[taxon.label].get('^ncbi:accession')
                if gi_id is None:
                    gi_id = self.otu_dict[taxon.label].get("^ot:originalLabel")
                new_label = "_".join([new_label, str(gi_id)])
                sp_counter = 2
                if new_label in new_names and norepeats:
                    # debug(self.otu_dict[taxon.label].keys())
                    # debug(gi_id)
                    new_label = "_".join([new_label, str(sp_counter)])
                    sp_counter += 1
            else:
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, taxon.label])
                    # debug(new_label)
            taxon.label = new_label
            new_names.add(new_label)
        tmp_tre.write(path=treepath,
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path=alnpath,
                      schema="fasta")

    def write_otus(self, filename, schema='table'):
        """Writes out OTU dict as json.

        :param filename: filename
        :param schema: either table or json format
        :return: writes out otu_dict to file
        """
        assert schema in ['table', 'json']
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln and tre and updates otu_dict,
        takes a single taxon_label as input.

        note: has test, test_remove_taxa_aln_tre.py

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return: removes information/data from taxon_label
        """
        debug('remove_taxa_aln_tre')
        # debug(taxon_label)
        # debug(type(taxon_label))
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        # debug(tax)
        # debug(tax2)
        # debug(len(self.tre.taxon_namespace))
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found
            # the first prune does not remove it sometimes...
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            # next line cannot happen here, as then the seq_dict_build_ crashes, if taxon was just added?
            # self.tre.taxon_namespace.remove_taxon_label(taxon_label)
            # debug(len(self.tre.taxon_namespace))
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, updated otu_dict entry but was never in tre or aln!"

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
    nexson = phy.get_study(study_id)['data']
    return nexson


def get_mrca_ott(ott_ids):
    """finds the mrca of the taxa in the ingroup of the original
    tree. The blast search later is limited to descendants of this
    mrca according to the ncbi taxonomy

    Used in the functions that generate the ATT object.

    :param ott_ids: list of all OToL identifiers for tiplabels in phylogeny
    :return: OToL identifier of most recent common ancestor or ott_ids
    """
    debug("get_mrca_ott")
    # drop_tip = []
    if None in ott_ids:
        ott_ids.remove(None)
    synth_tree_ott_ids = []
    ott_ids_not_in_synth = []
    for ott in ott_ids:
        # debug(ott)
        try:
            tree_of_life.mrca(ott_ids=[ott], wrap_response=False)
            synth_tree_ott_ids.append(ott)
        except:  # TODO: seems to be requests.exceptions.HTTPError: 500, don't know how to implement them
            debug("except")
            ott_ids_not_in_synth.append(ott)
            # drop_tip.append(ott)
    # debug(synth_tree_ott_ids)
    if len(synth_tree_ott_ids) == 0:
        sys.stderr.write('No sampled taxa were found in the current synthetic tree. '
                         'Please find and input and appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit()
    mrca_node = tree_of_life.mrca(ott_ids=synth_tree_ott_ids, wrap_response=False)  # need to fix wrap eventually
    # debug(mrca_node)
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
        sys.exit()
    return tax_id


def get_ott_ids_from_otu_dict(otu_dict):  # TODO put into data obj?
    """Get the ott ids from an otu dict object"""
    # TODO: never used
    ott_ids = []
    for otu in otu_dict:
        try:
            ott_ids.append(otu['^ot:ottId'])
        except KeyError:
            pass


#####################################

class IdDicts(object):
    """Wraps up the annoying conversions

    Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

    To build the class the following is needed:
        config_obj: Object of class config (see above)
        workdir: the path to the assigned working directory

    During the initializing process the following self objects are generated:
        self.workdir: contains path of working directory
        self.config: contains the Config class object
        self.ott_to_ncbi: dictionary
                        key: OToL taxon identifier
                        value: ncbi taxon identifier
        self.ncbi_to_ott: dictionary
                        key: OToL taxon identifier
                        value: ncbi taxon identifier
        self.ott_to_name: dictionary
                        key: OToL taxon identifier
                        value: OToL taxon name
        self.gi_ncbi_dict: dictionary
                        key: Genbank identifier
                        value: ncbi taxon identifier
        self.spn_to_ncbiid: dictionary
                        key: OToL taxon name
                        value: ncbi taxon identifier
        self.ncbiid_to_spn: dictionary
                        key: ncbi taxon identifier
                        value: ncbi taxon name
        Optional:
        self.ncbi_parser: initializes the ncbi_parser class, that contains information about rank and identifiers
        self.otu_rank: contains hierarchical information for web queries
    """

    # TODO - could - should be shared acrosss runs?! .... nooo.
    def __init__(self, config_obj, workdir, mrca = None):
        """Generates a series of name disambiguation dicts"""
        self.workdir = workdir  # TODO: Not needed. only used for dump and map_gi. map_gi file does not exists. dump is only used in wrapper, and we have the information of workdir available in wrapper functions anyways
        self.config = config_obj
        assert self.config.email
        self.ott_to_ncbi = {}  # currently only used to find mcra ncbi id from mrca ott_id
        self.ncbi_to_ott = {}  # used to get ott_id for new Genbank query taxa
        self.ott_to_name = {}  # used in add_otu to get name from otuId
        self.gi_ncbi_dict = {}  # file id_map doesn't exist (see below), is only filled by ncbi_parser (by subprocess in earlier versions of the code).
        self.spn_to_ncbiid = {}  # spn to ncbi_id, it's only fed by the ncbi_data_parser, but makes it faster
        self.ncbiid_to_spn = {}
        self.mrca_ott = mrca  # mrca_list
        self.mrca_ncbi = set()  # corresponding ids for mrca_ott list
        fi = open(config_obj.ott_ncbi)  # TODO need to keep updated, where does the file come from?
        for lin in fi:  # TODO This is insanely memory inefficient, how about using a pandas dataframe?
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()
            assert len(self.ott_to_ncbi) > 0
            assert len(self.ncbi_to_ott) > 0
            assert len(self.ott_to_name) > 0
        fi.close()
        # TODO: pandas solution? requires to rewrite usages of self.ott_to_ncbi, self.ncbi_to_ott, self.ott_to_name
        # assert os.path.exists(config_obj.ott_ncbi)
        # df = pd.read_csv(config_obj.ott_ncbi, sep='|', header=None, index_col=False,
        #                  names=[
        #                      'ott_id',
        #                      'ncbi_id',
        #                      'ott_name'
        #                  ])
        # assert len(df) > 0
        # TODO: where do we generate id_map?
        if os.path.isfile("{}/id_map.txt".format(workdir)):  # todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]
        if config_obj.blast_loc == 'remote':
            self.otu_rank = {}  # used only for web queries - contains taxonomic hierarchy information
        else:  # ncbi parser contains information about spn, tax_id, and ranks
            self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                                       nodes_file=self.config.ncbi_parser_nodes_fn)
        if self.mrca_ott != None:
            self.get_ncbi_mrca()

    def get_ncbi_mrca(self):
        """ get the ncbi tax ids from a list of mrca.
        """
        debug(type(self.mrca_ott))
        if type(self.mrca_ott) is not int:
            for ott_id in self.mrca_ott:
                debug(ott_id)
                debug(type(ott_id))
                self.ott_id_to_ncbiid(ott_id)
        else:
            self.ott_id_to_ncbiid(self.mrca_ott)
        debug(self.mrca_ncbi)

    def ott_id_to_ncbiid(self, ott_id):
        """ Find ncbi Id for ott id. Is only used for the mrca list thing.
        """
        if ott_id in self.ott_to_name:
            # debug("get id from ott_to_name")
            ott_name = self.ott_to_name[ott_id]
            if self.config.blast_loc == 'remote':
                debug(ott_name)
                self.get_rank_info_from_web(taxon_name=ott_name)
                debug(self.otu_rank.keys())
                ncbi_id = self.otu_rank[ott_name]["taxon id"]
            else:
                ncbi_id = self.ncbi_parser.get_id_from_name(ott_name)

        else:
            debug("else")
            tx = APIWrapper().taxomachine
            nms = tx.taxon(ott_id)
            debug(nms)
            ott_name = nms[u'unique_name']
            debug(ott_name)
            ncbi_id = None
            if u'ncbi' in nms[u'tax_sources']:
                ncbi_id = nms[u'tax_sources'][u'ncbi']
                debug(ncbi_id)
        if ncbi_id is not None:
            self.mrca_ncbi.add(ncbi_id)

    def get_ncbiid_from_tax_name(self, tax_name):
        """Get the ncbi_id from the species name using ncbi web query.

        :param tax_name: species name
        :return: corresponding ncbi id
        """
        ncbi_id = None
        if tax_name in self.spn_to_ncbiid:
            ncbi_id = self.spn_to_ncbiid[tax_name]
        else:
            
            debug(tax_name)
            try:
                debug("try2")
                tries = 15
                for i in range(tries):
                    # debug(Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0])
                    try:
                        Entrez.email = self.config.email
                        if tries >= 5:
                            ncbi_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
                        else:
                            tax_name = "'{}'".format(tax_name)
                            ncbi_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]

                        ncbi_id = int(ncbi_id)
                    except:  # TODO: is either IndexError or rllib2.HTTPError: HTTP Error 400: Bad Request
                        # debug("except esearch/read")
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            raise
                    break
                # debug(ncbi_id)
                # debug(type(ncbi_id))
            except:  # TODO: is either IndexError or rllib2.HTTPError: HTTP Error 400: Bad Request
                debug("except")
                try:
                    ncbi = NCBITaxa()
                    tax_info = ncbi.get_name_translator([tax_name])
                    debug(tax_info)
                    if tax_info == {}:
                        tax_name = "'{}'".format(tax_name)
                        tax_info = ncbi.get_name_translator([tax_name])
                    ncbi_id = int(tax_info.items()[0][1][0])
          	except: 
                    sys.stderr.write("Taxon name does not match any name in ncbi. Check that name is written "
                                     "correctly: {}! We set it to unidentified".format(tax_name))
                    tax_name = 'Eukaryota'
                    ncbi_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
                    tax_name = 'unidentified_Eukaryota'
                    ncbi_id = int(ncbi_id)
        assert type(ncbi_id) is int
        self.spn_to_ncbiid[tax_name] = ncbi_id
        return ncbi_id

    def get_rank_info_from_web(self, taxon_name):
        # NOTE MK: old version for web-queries, needs to be implemented newly
        """Collects rank and lineage information from ncbi,
        used to delimit the sequences from blast,
        when you have a local blast database or a Filter Blast run
        """
        # debug("get_rank_info")
        # Entrez.email = self.config.email
        # if gi_id:
        #     # debug("gi_id to tax_name")
        #     tries = 5
        #     for i in range(tries):
        #         try:
        #             handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
        #         except:
        #             if i < tries - 1:  # i is zero indexed
        #                 continue
        #             else:
        #                 raise
        #         break
        #     read_handle = Entrez.read(handle)[0]
        #     tax_name = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value']
        # else:
        #     tax_name = str(taxon_name).replace("_", " ")
        tax_name = taxon_name.replace(" ", "_")
        if tax_name not in self.otu_rank.keys():

            ncbi_id = self.get_ncbiid_from_tax_name(tax_name)
            # replaced by newer functio: get_ncbiid_from_tax_name()
            # debug("tax_name to rank")
            # ncbi = NCBITaxa()
            # try:
            #     tax_id = int(Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0])
            # except:
            #     # debug("except")
            #     tax_info = ncbi.get_name_translator([tax_name])
            #     if tax_info == {}:
            #         print("Taxon name does not match any species name in ncbi. Check that name is written correctly!")
            #     tax_id = int(tax_info.items()[0][1][0])
            ncbi = NCBITaxa()
            lineage = ncbi.get_lineage(ncbi_id)
            lineage2ranks = ncbi.get_rank(lineage)
            tax_name = str(tax_name).replace(" ", "_")
            assert type(ncbi_id) == int
            self.otu_rank[tax_name] = {"taxon id": ncbi_id, "lineage": lineage, "rank": lineage2ranks}
        return tax_name


    def find_name(self, sp_dict=None, gi=None):
        """ Find the taxon name in the sp_dict or of a gi_id.
        If not already known it will ask ncbi using the gi_id.
        """
        debug("find_name")
        inputinfo = False
        if sp_dict is not None or gi is not None:
            inputinfo = True
        assert inputinfo is True
        tax_name = None
        # debug([sp_dict, gi])
        if sp_dict:
            # debug('^ot:ottTaxonName' in sp_dict)
            # debug(u'^ot:ottTaxonName' in sp_dict)
            # debug(sp_dict.keys())
            if '^ot:ottTaxonName' in sp_dict:
                tax_name = sp_dict['^ot:ottTaxonName']
            elif '^user:TaxonName' in sp_dict:
                tax_name = sp_dict['^user:TaxonName']
        if tax_name is None:
            gi_id = None
            if gi:
                gi_id = gi
            elif '^ncbi:gi' in sp_dict:
                gi_id = sp_dict['^ncbi:gi']
                if gi_id is None:
                    gi_id = sp_dict['^ncbi:accession']

            else:
                sys.stderr.write("There is no name supplied and no gi available. This should not happen! Check name!")
            if gi_id in self.gi_ncbi_dict:
                ncbi_id = self.gi_ncbi_dict[gi_id]
                # # deprecated for ncbi_parser
                # for key, value in self.spn_to_ncbiid.values():
                #     if value == ncbi_id:
                #         tax_name = key
                if ncbi_id in self.ncbiid_to_spn.keys():
                    tax_name = self.ncbiid_to_spn[ncbi_id]
                else:
                    tax_name = self.ncbi_parser.get_name_from_id(ncbi_id)
                    self.ncbiid_to_spn[ncbi_id] = tax_name
            else:  # usually being used for web-queries, local blast searches should have the information
                # debug(gi_id)
                # debug(type(gi_id))
                tries = 10
                Entrez.email = self.config.email
                for i in range(tries):
                    try:
                        # debug("find name efetch")
                        handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
                    except IndexError:
                        # debug("except efetch")
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            raise
                    break
                read_handle = Entrez.read(handle)
                handle.close()
                tax_name = get_ncbi_tax_name(read_handle)
                ncbi_id = get_ncbi_tax_id(read_handle)
                self.ncbiid_to_spn[ncbi_id] = tax_name
                self.gi_ncbi_dict[gi_id] = ncbi_id
                if sp_dict:
                    sp_dict['^ot:ottTaxonName'] = tax_name
                    sp_dict['^ncbi:taxon'] = ncbi_id
        assert tax_name is not None
        tax_name = tax_name.replace(" ", "_")
        return tax_name

    def map_gi_ncbi(self, gi_id):
        """get the ncbi taxon id's for a gi input.

        Finds different identifiers and information from a given gi_id and fills the corresponding self.objects
        with the retrieved information.

        :param gi_id: Genbank ID to
        :return: ncbi taxon id
        """
        # debug("map_gi_ncbi")
        if _DEBUG == 2:
            sys.stderr.write("mapping gi {}\n".format(gi_id))
        if gi_id in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi_id])
        else:
            # debug(gi)
            tax_name = self.find_name(gi=gi_id)
            if self.config.blast_loc == 'remote':
                try:
                    self.get_rank_info_from_web(taxon_name=tax_name)
                    tax_id = self.otu_rank[tax_name]["taxon id"]
                except IndexError:  # get id via genbank query xref in description
                    tries = 10
                    Entrez.email = self.config.email
                    for i in range(tries):
                        try:
                            # debug("find name efetch")
                            handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
                        except IndexError:
                            # debug("except efetch")
                            if i < tries - 1:  # i is zero indexed
                                continue
                            else:
                                raise
                        break
                    read_handle = Entrez.read(handle)
                    handle.close()
                    tax_name = get_ncbi_tax_name(read_handle)
                    ncbi_id = get_ncbi_tax_id(read_handle)
            else:
                tax_id = self.ncbi_parser.get_id_from_name(tax_name)
            self.ncbiid_to_spn[tax_id] = tax_name
            self.gi_ncbi_dict[gi_id] = tax_id
            self.spn_to_ncbiid[tax_name] = tax_id
        return tax_id

    def dump(self, filename=None):
        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/id_pickle.p".format(self.workdir, filename), "wb")
        pickle.dump(self, ofi)


class PhyscraperScrape(object):  # TODO do I want to be able to instantiate this in a different way?!
    # set up needed variables as nones here?!
    # TODO better enforce ordering
    """This is the class that does the perpetual updating

        To build the class the following is needed:
            data_obj: Object of class ATT (see above)
            ids_obj: Object of class IdDict (see above)

        During the initializing process the following self.objects are generated:
            self.workdir: path to working directory retrieved from ATT object = data_obj.workdir
            self.logfile: path of logfile
            self.data = ATT object
            self.ids = IdDict object
            self.config = Config object
            * value: is a dictionary again, but varies depending on the blast method used:
                * local blast:
                    * ^ncbi:gi': gi id
                    * 'accession':Genbank accession number
                    * 'staxids': ncbi taxon id
                    * 'sscinames': ncbi species name
                    * 'pident': pident
                    * 'evalue': evalue
                    * 'bitscore': bitscore
                    * 'sseq': sequence
                    * 'title': title from GenBank sequence
                * web-query blast:
                    * 'accession':Genbank accession number
                    * 'length': length of sequence
                    * 'title': string combination of hit_id and hit_def
                    * 'hit_id': string combination of gi id and accession number
                    * 'hsps': Bio.Blast.Record.HSP object
                    * 'hit_def': title from GenBank sequence
            self.otu_by_gi: dictionary that contains ????:
                        key:
                        value:
            self._to_be_pruned: list that contains ????
            self.mrca_ncbi: ncbi identifier of mrca

            self.tmpfi: path to a file or folder???
            self.blast_subdir: path to folder that contains the files writen during blast

            self.newseqs_file = filename of files that contains the sequences from self.new_seqs_otu_id
            self.date: Date of the run - may lag behind real date!
            self.repeat: either 1 or 0, it is used to determine if we continue updating the tree, no new seqs found = 0
            self.newseqsgi: List of all gi_ids that were passed to remove_identical_seq().
                            It is used to speed up adding process
            self.blacklist: list of gi_id of sequences that shall not be added or need to be removed. Supplied by user.
            self.gi_list_mrca: List of all gi_ids available on GenBank for a given mrca.
                               Used to limit possible seq to add.
            self.seq_filter: List of words that may occur in otu_dict.status and which shall not be used for
                            building the FilterBlast.sp_d (that's the main function), but it is also used as assert
                            statement to make sure unwanted seqs are not added.
            self.unpublished: True/False. Used to look for local unpublished seq that shall be added if True.
            self.path_to_local_seq: Usually False, contains path to unpublished sequences if option is used.

        Following functions are called during the init-process:
            self.reset_markers():
                adds things to self:  # TODO: I think they are used to make sure certain function run,
                                     if program crashed and pickle file if read in.
                    self._blasted: 0/1, if run_blast() was called, it is set to 1 for the round.
                    self._blast_read: 0/1, if read_blast() was called, it is set to 1 for the round.
                    self._identical_removed: 0
                    self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
                    self._query_seqs_aligned: 0
                    self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
                    self._reconciled: 0
                    self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
    """

    def __init__(self, data_obj, ids_obj):
        # todo check input types assert()
        debug("start base class init")
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = self.ids.config  # TODO: this is already part of self.ids, information are doubled.
        self.new_seqs = {}  # all new seq after read_blast
        self.new_seqs_otu_id = {}  # only new seq which passed remove_identical
        self.otu_by_gi = {}  # TODO: What was this intended for?
        self._to_be_pruned = []  # TODO: What was this intended for? We don't use it
        self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir) # TODO: For what do we want to use this?
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "tmp.fasta"  # TODO: is renamed before using it. Name is easily confused with 'tmp.fas' which is used as temporary file that contains the sequence that is currently blasted
        self.date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        self.repeat = 1  # used to determine if we continue updating the tree
        self.newseqsgi = []  # all ever added gi during any PhyScraper run, used to speed up adding process
        self.blacklist = []  # remove sequences by default
        self.gi_list_mrca = []  # all gi_ids of a given mrca. Used to limit possible seq to add. 
        if self.config.blast_loc == 'local' and len(self.gi_list_mrca) == 0:
           self.gi_list_mrca = self.get_all_gi_mrca()
            # debug(self.gi_list_mrca)
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,", "local"]  # TODO MK: try to move completely to FilterBlast class
        self.reset_markers()
        self.unpublished = False  # used to look for local unpublished seq that shall be added.
        self.path_to_local_seq = False  # path to unpublished seq.
        self.backbone = False
        self.OToL_unmapped_tips()  # added to remove un-mapped tips from OToL
        self.ids.ingroup_mrca = data_obj.ott_mrca

###############################3


        if _deep_debug == 1:
            self.newadd_gi_otu = {}  # search for doubles!

    # TODO is this the right place for this?
    def reset_markers(self):
        self._blasted = 0
        self._blast_read = 0
        self._identical_removed = 0  # TODO: We don't use it
        self._query_seqs_written = 0
        self._query_seqs_aligned = 0  # TODO: We don't use it
        self._query_seqs_placed = 0
        self._reconciled = 0  # TODO: We don't use it
        self._full_tree_est = 0

    def OToL_unmapped_tips(self):
        """Assign names or remove tips from aln and tre that were not mapped during initiation of ATT class.
        """
        debug("OTOL unmapped")
        # debug(self.config.unmapped)
        # debug(self.data.ott_mrca)
        # debug(len(self.data.aln.taxon_namespace))
        # debug(self.data.aln.taxon_namespace)
        if self.config.unmapped == 'remove':
            # debug("remove_OToL_unmapped")
            # drop tips without ott _id
            # debug(ott_ids_not_in_synth)
            # debug(len(self.data.otu_dict))
            # debug(len(self.aln.taxon_namespace))
            for key in self.data.otu_dict:
                # debug(self.data.otu_dict[key].keys())
                if '^ot:ottId' not in self.data.otu_dict[key]:
                    # second condition for OToL unmapped taxa, not present in own_data
                    if u'^ot:treebaseOTUId' in self.data.otu_dict[key]:
                        # debug(key)
                        self.data.remove_taxa_aln_tre(key)
            # debug(len(self.data.otu_dict))
            # # debug(self.aln.taxon_namespace)        
            # debug(len(self.data.aln.taxon_namespace))
            # # debug(self.tre.taxon_namespace)        
            # debug(len(self.data.tre.taxon_namespace))
        else:
            # debug("keep_OToL_unmapped")
            i=1
            for key in self.data.otu_dict:
                i=i+1
                # debug(i)
                # debug(self.data.otu_dict[key].keys())
                if '^ot:ottId' not in self.data.otu_dict[key]:
                    # debug("no id known")
                    self.data.otu_dict[key]['^ot:ottId'] = self.data.ott_mrca
                    if self.data.ott_mrca in self.ids.ott_to_name:
                        # debug("get id from ott_to_name")
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
        cwd = os.getcwd()
        os.chdir(self.config.blastdb)
        toblast = open("{}/tmp.fas".format(self.blast_subdir), 'w')
        toblast.write(">{}\n".format(taxon_label))
        toblast.write("{}\n".format(query))
        toblast.close()
        # this formats allows to get the taxonomic information at the same time
        outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq stitle'"
        # outfmt = " -outfmt 5"  # format for xml file type
        # TODO query via stdin
        # TODO MK: update to blast+ v. 2.8 - then we can limit search to taxids: -taxids self.mrca_ncbi
        blastcmd = "blastn -query " + \
                   "{}/tmp.fas".format(self.blast_subdir) + \
                   " -db {}nt -out ".format(self.config.blastdb) + \
                   fn_path + \
                   " {} -num_threads {}".format(outfmt, self.config.num_threads) + \
                   " -max_target_seqs {} -max_hsps {}".format(self.config.hitlist_size,
                                                              self.config.hitlist_size)
        # debug(blastcmd)
        os.system(blastcmd)
        os.chdir(cwd)

    def run_web_blast_query(self, query, equery, fn_path):
        """Equivalent to run_local_blast_cmd() but for webqueries, 
        that need to be implemented differently.

        :param query: query sequence
        :param equery: method to limit blast query to mrca
        :param fn_path: path to output file for blast query result
        :return: runs web blast query and writes it to file
        """
        if self.config.url_base:
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          url_base=self.config.url_base,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size,
                                          num_threads=self.config.num_threads)
        else:
            debug("use BLAST webservice")
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size)
            # debug(result_handle.read())
        save_file = open(fn_path, "w")
        save_file.write(result_handle.read())
        result_handle.close()
        save_file.close()

    def run_blast(self, delay=14):  # TODO Should this be happening elsewhere?
        """generates the blast queries and saves them depending on the blasting method to different file formats

        :param delay: number that determines when a previously blasted sequence is reblasted - time is in days
        :return: writes blast queries to file
        """
        debug("run_blast")
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            # debug("{}, {}".format(taxon, seq))
            # TODO temp until I fix delete
            if otu_id in self.data.otu_dict:
                if _VERBOSE:
                    sys.stdout.write("blasting {}\n".format(otu_id))
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days)
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                if self.unpublished:
                    debug("run against local unpublished data")
                    toblast = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                    toblast.write(">{}\n".format(taxon.label))
                    toblast.write("{}\n".format(query))
                    toblast.close()
                    blast_db = "local_unpubl_seq_db"
                    output = "tst_fn"
                    blastcmd = "blastn -query {}/tmp.fas -db {} -out output_{}.xml " \
                               "-outfmt 5".format(self.blast_subdir, blast_db, output)
                    os.system(blastcmd)
                    if self.backbone == True:
                        self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                    # self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                else:
                    if time_passed > delay:
                        if self.config.blast_loc == 'local':
                            file_ending = "txt"
                        else:
                            file_ending = "xml"
                        if self.config.gifilename is True:
                            fn = self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label)
                            if fn == None:
                                fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label)

                            fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                        else:
                            fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                        if _DEBUG:
                            sys.stdout.write("attempting to write {}\n".format(fn_path))
                        if not os.path.isfile(fn_path):
                            if _VERBOSE:
                                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                            if self.config.blast_loc == 'local':
                                self.run_local_blast_cmd(query, taxon.label, fn_path)
                                # cwd = os.getcwd()
                                # os.chdir(self.config.blastdb)
                                # toblast = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                                # toblast.write(">{}\n".format(taxon.label))
                                # toblast.write("{}\n".format(query))
                                # toblast.close()
                                # # this formats allows to get the taxonomic information at the same time
                                # outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq stitle'"
                                # # outfmt = " -outfmt 5"  # format for xml file type
                                # # TODO query via stdin
                                # blastcmd = "blastn -query " + \
                                #            "{}/tmp.fas".format(self.blast_subdir) + \
                                #            " -db {}nt -out ".format(self.config.blastdb) + \
                                #            fn_path + \
                                #            " {} -num_threads {}".format(outfmt, self.config.num_threads) + \
                                #            " -max_target_seqs {} -max_hsps {}".format(self.config.hitlist_size,
                                #                                                       self.config.hitlist_size)
                                # # debug(blastcmd)
                                # os.system(blastcmd)
                                # os.chdir(cwd)
                                # self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                            if self.config.blast_loc == 'remote':
                                debug(len(self.ids.mrca_ncbi))

                                if len(self.ids.mrca_ncbi) >=2:
                                    len_ncbi = len(self.ids.mrca_ncbi)
                                    equery = ''
                                    for ncbi_id in self.ids.mrca_ncbi:
                                        if len_ncbi >= 2:
                                            equery = equery + "txid{}[orgn] OR ".format(ncbi_id)
                                            len_ncbi = len_ncbi - 1
                                        else:
                                            equery = equery + "txid{}[orgn]) ".format(ncbi_id)
                                    debug(equery)
                                    equery = "(" + equery + "AND {}:{}[mdat]".format(last_blast, today)
                                    debug(equery)
                                else:
                                    equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi, last_blast, today)
                                    self.run_web_blast_query(query, equery, fn_path)
                                # if self.config.url_base:
                                #     result_handle = AWSWWW.qblast("blastn",
                                #                                   "nt",
                                #                                   query,
                                #                                   url_base=self.config.url_base,
                                #                                   entrez_query=equery,
                                #                                   hitlist_size=self.config.hitlist_size,
                                #                                   num_threads=self.config.num_threads)
                                # else:
                                #     debug("use BLAST webservice")
                                #     result_handle = AWSWWW.qblast("blastn",
                                #                                   "nt",
                                #                                   query,
                                #                                   entrez_query=equery,
                                #                                   hitlist_size=self.config.hitlist_size)
                                #     # debug(result_handle.read())
                                # save_file = open(fn_path, "w")
                                # save_file.write(result_handle.read())
                                # result_handle.close()
                                # save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                        #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
                        else:
                            # changes date of blasted accordingly, if file is already present in the folder
                            if _DEBUG:
                                sys.stdout.write("file {} exists in current blast run. Will not blast, "
                                                 "delete file to force\n".format(fn_path))
                            if _DEBUG_MK == 1:
                                self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                    else:
                        if _VERBOSE:
                            sys.stdout.write("otu {} was last blasted {} days ago and is not being re-blasted. "
                                             "Use run_blast(delay = 0) to force a search.\n".format(otu_id, last_blast))
        self._blasted = 1

    def get_all_gi_mrca(self):
        """get all available gi numbers from Genbank for mrca.

        The list will be used to filter out sequences from the local Blast search,
        that do not belong to ingroup

        :return: list of corresponding gi numbers
        """
        # TODO MK: gi list limited to 100000000, for huge trees that is a problem. Think of what to do...
        debug("get_all_gi_mrca")
        Entrez.email = self.config.email
### NOTE MK> if mrca ingroup is list, use the list not the mrca_ncbi/ott
        if len(self.ids.mrca_ncbi) >=2:
            len_ncbi = len(self.ids.mrca_ncbi)
            equery = '('
            for ncbi_id in self.ids.mrca_ncbi:
                if len_ncbi >= 2:
                    equery = equery + "txid{}[orgn] OR ".format(ncbi_id)
                    len_ncbi = len_ncbi - 1
                else:
                    equery = equery + "txid{}[orgn])".format(ncbi_id)
            debug(equery)
        else:
            equery = "txid{}[orgn]".format(self.mrca_ncbi)
        handle = Entrez.esearch(db="nucleotide", term=equery,
                                usehistory='n', RetMax=100000000)
        records = Entrez.read(handle)
        id_list = records['IdList']
        id_list = [int(x) for x in id_list]
        return id_list

    def read_local_blast(self, fn_path):
        """ Implementation to read in results of local blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gi_dict dictionaries
        """
        query_dict = {}
        with open(fn_path, mode='r') as infile:
            for lin in infile:
                # print(lin)
                sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, stitle = lin.strip().split('\t')
                gi_id = int(sseqid.split("|")[1])
                gi_acc = sseqid.split("|")[3]
                sseq = sseq.replace("-", "")
                # sscinames = sscinames.replace(" ", "_").replace("/", "_").replace("-", "_")
                sscinames = sscinames.replace(" ", "_").replace("/", "_")
                # debug(staxids)
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                # print(pident, evalue, bitscore)
                # print(type(pident), type(evalue), type(bitscore))
                # NOTE: sometimes there are seq which are identical and are combined in the local blast db, just get first one
                if len(staxids.split(";")) > 1:
                    staxids = int(staxids.split(";")[0])
                    # debug(staxids)
                    sscinames = sscinames.split(";")[0]
                else:
                    staxids = int(staxids)
                    # debug(sscinames)
                # print(type(staxids))
                assert type(staxids) is int
                self.ids.spn_to_ncbiid[sscinames] = staxids
                if gi_id not in self.ids.gi_ncbi_dict:  # fill up dict with more information.
                    self.ids.gi_ncbi_dict[gi_id] = staxids
                if gi_id not in query_dict and gi_id not in self.newseqsgi:
                    query_dict[gi_id] = {'^ncbi:gi': gi_id, 'accession': gi_acc, 'staxids': staxids,
                                         'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                         'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
        for key in query_dict.keys():
            if float(query_dict[key]['evalue']) < float(self.config.e_value_thresh):
                gi_id = query_dict[key]['^ncbi:gi']
                if gi_id == None:
                    gi_id = query_dict[key]['accession']

                # debug(type(gi_id))
                # debug(gi_id)
                if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                    # debug("pass")
                    pass
                else:
                    # debug("try to add to new seqs")
                    if gi_id not in self.data.gi_dict:  # skip ones we already have            
                        # debug("added")
                        self.new_seqs[gi_id] = query_dict[key]['sseq']
                        self.data.gi_dict[gi_id] = query_dict[key]

    def read_webbased_blast_query(self, fn_path):
        """ Implementation to read in results of web blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gi_dict dictionaries
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
                            gi_id = int(alignment.title.split('|')[1])
                            assert type(gi_id) is int
                            if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                                pass
                            else:
                                if gi_id not in self.data.gi_dict:  # skip ones we already have
                                    # debug("add gi to new seqs")
                                    self.new_seqs[gi_id] = hsp.sbjct
                                    self.data.gi_dict[gi_id] = alignment.__dict__
        except ValueError:
            sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

    def read_blast(self, blast_dir=None):
        """reads in and processes the blast xml files

        :param blast_dir: path to directory which contains blast files
        :return: fills different dictionaries with information from blast files
        """
        debug("read blast")
        # debug(blast_dir)
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
            # use read_local_blast func
            output_blast = "output_tst_fn.xml"
            gi_counter = 1
            general_wd = os.getcwd()
            os.chdir(os.path.join(self.workdir, "blast"))
            xml_file = open(output_blast)
            os.chdir(general_wd)
            blast_out = NCBIXML.parse(xml_file)
            for blast_record in blast_out:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if float(hsp.expect) < float(self.config.e_value_thresh):
                            gi_id = alignment.title.split('|')[-1].split(' ')[-1]
                            # debug(gi_id)
                            if gi_id not in self.data.gi_dict:  # skip ones we already have
                                # debug("add gi to new seqs")
                                self.make_otu_dict_entry_unpubl(gi_id)
                                fake_gi = "unpubl_{}".format(gi_id)
                                self.new_seqs[fake_gi] = hsp.sbjct
                                # debug(gi_id)
                                # print(self.data.unpubl_otu_json)
                                # print(self.data.unpubl_otu_json.keys())

                                # print(self.data.unpubl_otu_json['otu{}'.format(gi_id)])
                                self.data.gi_dict[fake_gi] = {'accession': "000000{}".format(gi_counter),
                                                              'title': "unpublished", 'localID': gi_id}
                                self.data.gi_dict[fake_gi].update(self.data.unpubl_otu_json['otu{}'.format(gi_id)])
                                gi_counter += 1
                                # self.data.gi_dict[fake_gi] = alignment.__dict__
        else:
            if not self._blasted:
                self.run_blast()
            assert os.path.exists(self.blast_subdir)
            for taxon in self.data.aln:
                # debug(taxon)
                # debug("add blast seq to new seqs")
                # debug("blast location is: {}".format(self.config.blast_loc))
                if self.config.blast_loc == 'local':
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                # debug(self.config.gifilename)
                # debug(self.blast_subdir)
                if self.config.gifilename is True:
                    fn = self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label)
                    if fn == None:
                        fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label)

                    fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                else:
                    fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                if _DEBUG:
                    sys.stdout.write("attempting to read {}\n".format(fn_path))
                if os.path.isfile(fn_path):
                    if self.config.blast_loc == 'local':  # new method to read in txt format
                        self.read_local_blast(fn_path)
                        # moved to separate function
                        # query_dict = {}
                        # with open(fn_path, mode='r') as infile:
                        #     for lin in infile:
                        #         # print(lin)
                        #         sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, stitle = lin.strip().split('\t')
                        #         gi_id = int(sseqid.split("|")[1])
                        #         gi_acc = sseqid.split("|")[3]
                        #         sseq = sseq.replace("-", "")
                        #         # sscinames = sscinames.replace(" ", "_").replace("/", "_").replace("-", "_")
                        #         sscinames = sscinames.replace(" ", "_").replace("/", "_")
                                
                        #         # debug(staxids)
                        #         pident = float(pident)
                        #         evalue = float(evalue)
                        #         bitscore = float(bitscore)
                        #         # print(pident, evalue, bitscore)
                        #         # print(type(pident), type(evalue), type(bitscore))
                        #         # NOTE: sometimes there are seq which are identical and are combined in the local blast db, just get first one
                        #         if len(staxids.split(";")) > 1:
                        #             staxids = int(staxids.split(";")[0])
                        #             # debug(staxids)
                        #             sscinames = sscinames.split(";")[0]
                        #         else:
                        #             staxids = int(staxids)
                        #             # debug(sscinames)
                        #         # print(type(staxids))
                        #         assert type(staxids) is int
                        #         self.ids.spn_to_ncbiid[sscinames] = staxids
                        #         if gi_id not in self.ids.gi_ncbi_dict:  # fill up dict with more information.
                        #             self.ids.gi_ncbi_dict[gi_id] = staxids
                        #         if gi_id not in query_dict and gi_id not in self.newseqsgi:
                        #             query_dict[gi_id] = {'^ncbi:gi': gi_id, 'accession': gi_acc, 'staxids': staxids,
                        #                                  'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                        #                                  'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                        # for key in query_dict.keys():
                        #     if float(query_dict[key]['evalue']) < float(self.config.e_value_thresh):
                        #         gi_id = query_dict[key]['^ncbi:gi']
                        #         # debug(type(gi_id))
                        #         # debug(gi_id)
                        #         if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                        #             # debug("pass")
                        #             pass
                        #         else:
                        #             # debug("try to add to new seqs")
                        #             if gi_id not in self.data.gi_dict:  # skip ones we already have            
                        #                 # debug("added")
                        #                 self.new_seqs[gi_id] = query_dict[key]['sseq']
                        #                 self.data.gi_dict[gi_id] = query_dict[key]
                    else:
                        self.read_webbased_blast_query(fn_path)
                        # moved to own function.
                        # result_handle = open(fn_path)
                        # try:
                        #     if _VERBOSE:
                        #         sys.stdout.write(".")
                        #     blast_records = NCBIXML.parse(result_handle)
                        #     for blast_record in blast_records:
                        #         for alignment in blast_record.alignments:
                        #             for hsp in alignment.hsps:
                        #                 if float(hsp.expect) < float(self.config.e_value_thresh):
                        #                     gi_id = int(alignment.title.split('|')[1])
                        #                     assert type(gi_id) is int
                        #                     if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                        #                         pass
                        #                     else:
                        #                         if gi_id not in self.data.gi_dict:  # skip ones we already have
                        #                             # debug("add gi to new seqs")
                        #                             self.new_seqs[gi_id] = hsp.sbjct
                        #                             self.data.gi_dict[gi_id] = alignment.__dict__
                        # except ValueError:
                        #     sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

        self.date = str(datetime.date.today())
        debug("len new seqs dict after evalue filter")
        debug(len(self.new_seqs))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from GenBank after evalue filtering\n".format(len(self.new_seqs)))

        self._blast_read = 1

    # TODO this should go back in the class and should prune the tree

    def get_sp_id_of_otulabel(self, label):
        """Get the species name and the corresponding ncbi id of the otu.

        :param label: otu_label = key from otu_dict
        :return: ncbi id of corresponding label
        """
        debug("get_tax_id_of_otulabel")
        debug(label)
        # debug(self.data.otu_dict[label])
        spn_of_label = self.ids.find_name(sp_dict=self.data.otu_dict[label])
        debug(spn_of_label)
        if spn_of_label is not None:
            # spn_of_label = str(spn_of_label).replace(" ", "_").replace("-", "_")
            spn_of_label = str(spn_of_label).replace(" ", "_")
        else:
            debug("Problem, no tax_name found!")
        debug(self.data.otu_dict[label])
        if '^ncbi:taxon' in self.data.otu_dict[label]:
            id_of_label = self.data.otu_dict[label]['^ncbi:taxon']
        elif spn_of_label in self.ids.spn_to_ncbiid:
            id_of_label = self.ids.spn_to_ncbiid[spn_of_label]
        elif u'^ot:ottId' in self.data.otu_dict[label]:  # from OTT to ncbi id
            info = get_ott_taxon_info(spn_of_label.replace("_", " "))
            debug(info)
            ottid, ottname, id_of_label = info
            assert ottid == self.data.otu_dict[label][u'^ot:ottId']
            self.ids.ncbi_to_ott[id_of_label] = ottid
            self.ids.ott_to_name[ottid] = spn_of_label
        else:
            if self.config.blast_loc == 'remote':
                self.ids.get_rank_info_from_web(taxon_name=spn_of_label)
                id_of_label = self.ids.otu_rank[spn_of_label]["taxon id"]
            else:
                id_of_label = self.ids.ncbi_parser.get_id_from_name(spn_of_label)
            self.ids.spn_to_ncbiid[spn_of_label] = id_of_label
        return id_of_label

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
        # TODO unify spp name somehow?
        # debug("seq_dict_build")
        # debug(label)
        id_of_label = self.get_sp_id_of_otulabel(label)
        # debug("id of label is{}".format(id_of_label))
        
        ##############################################
        if _deep_debug == 1:
            label_gi_id = self.data.otu_dict[label]['^ncbi:gi']
            deep_debug("label_gi_id is{}".format(label_gi_id))
            for item in self.data.otu_dict:
                deep_debug("double-check")
                if '^ncbi:gi' in self.data.otu_dict[item]:
                    deep_debug("cehck doubles")
                    if self.data.otu_dict[item]['^ncbi:gi'] == label_gi_id and label != item:
                        exit(-1)
        
        #################################################
        # debug(id_of_label)
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        continue_search = False
        never_add = False
        for tax_lab in tax_list:
            existing_id = self.get_sp_id_of_otulabel(tax_lab)
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            # debug("length")
            # debug("{}, {}".format(len(inc_seq), len(new_seq)))
            # debug(sum(self.data.orig_seqlen) / len(self.data.orig_seqlen))
            # debug("{}, {}".format(id_of_label, existing_id))
            if len(new_seq) >= sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * 2.5:
                debug("seq not added because it's to long...")
            elif len(inc_seq) >= len(new_seq):  # if seq is identical and shorter
                if inc_seq.find(new_seq) != -1:
                    # if (existing_taxa != spn_of_label and existing_taxa is not None) or
                    # print(type(existing_id))
                    if type(existing_id) == int and existing_id != id_of_label:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, "
                                             "but different species name\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; " \
                                                                          "subsequence, but different species"
                        seq_dict[label] = seq
                        debug("{} and {} are subsequences, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:  # subseq of same sp.
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        debug("{} not added, subseq of {}".format(id_of_label, existing_id))
                        never_add = True
                        continue
                    # return seq_dict
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
                        # elif spn_of_label not in exists:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, but different "
                                             "species concept\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; supersequence, " \
                                                                          "but different species"
                        seq_dict[label] = seq
                        debug("{} and  {} supersequence, but different sp. concept".format(id_of_label, existing_id))
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
                        debug("{} added, instead of  {}".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    # return seq_dict

        if continue_search is True or never_add is True:
            if (self.data.otu_dict[label]['^physcraper:status'].split(' ')[0] in self.seq_filter) or never_add is True:
                if label in seq_dict.keys():
                    del seq_dict[label]
                # else:
                #     debug("label was never added to seq_dict")
                # try:
                if label in self.data.aln.taxon_namespace or label in self.data.tre.taxon_namespace:
                    self.data.remove_taxa_aln_tre(label)
                # except:
                else:
                    debug("label was never added to aln or tre")
                # Note: should not be the word 'deleted', as this is used in self.seq_filter
                self.data.otu_dict[label]['^physcraper:status'] = "removed in seq dict build"
                return seq_dict
        if _VERBOSE:
            sys.stdout.write(".")
            if i % 50 == 0:
                sys.stdout.write("\n")
        seq_dict[label] = seq

        #####################################
        if _deep_debug == 1:
            for item in self.data.otu_dict:
                # debug(self.data.otu_dict[item])
                if '^ncbi:gi' in self.data.otu_dict[item]:
                    debug("check doubles. {} vs old:".format(label_gi_id, self.data.otu_dict[item]['^ncbi:gi']))
                    debug("{}, {}".format(label, item))
                    # debug(some)

                    if self.data.otu_dict[item]['^ncbi:gi'] == label_gi_id and label != item:
                        exit(-1)
        ########################################
        return seq_dict

    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        """
        debug("remove identical seqs")
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        # Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        assert self.config.seq_len_perc <= 1
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        # debug(self.config.seq_len_perc)
        # debug("seq cutoff")
        # debug(seq_len_cutoff)
        for gi_id, seq in self.new_seqs.items():
            # debug(gi_id)
            if self.blacklist is not None and gi_id in self.blacklist:
                debug("gi_id in blacklist, not added")
                pass
            elif gi_id in self.newseqsgi:  # added to increase speed. often seq was found in another blast file
                debug("passed, was already added")
                pass
            else:
                # debug("add to aln if no similar seq exists")
                if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                    if type(gi_id) == int or gi_id.isdigit():
                        # debug("gi_id is digit")
                        if type(gi_id) != int:
                            sys.stdout.write("WARNING: gi_id {} is no integer. "
                                             "Will convert value to int\n".format(gi_id))
                            debug("WARNING: gi_id {} is no integer. Will convert value to int\n".format(gi_id))
                            gi_id = int(gi_id)
                            #########################
                            if _deep_debug == 1: 
                                currentgilist = self.find_otudict_gi()
                                debug(currentgilist)
                                if gi_id in currentgilist:
                                    exit(-1)
                                    ####################
                    self.newseqsgi.append(gi_id)
                    otu_id = self.data.add_otu(gi_id, self.ids)
                    # debug(otu_id)
                    # debug("go to seq_dict_build")
                    # debug(self.data.otu_dict[otu_id])

                    self.seq_dict_build(seq, otu_id, tmp_dict)
            # else:
            #     debug("gi was already compared")
        old_seqs_ids = set()
        for tax in old_seqs:
            old_seqs_ids.add(tax)
        assert old_seqs_ids.issubset(tmp_dict.keys())
        for tax in old_seqs:
           # try:
            del tmp_dict[tax]
           # except KeyError:
           #     pass
        self.new_seqs_otu_id = tmp_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
       
        ##################################
        if _deep_debug == 1:
            for otu in self.new_seqs_otu_id:
                otu_gi_id = self.data.otu_dict[otu]['^ncbi:gi']
                if otu_gi_id not in self.newadd_gi_otu:
                    self.newadd_gi_otu[otu_gi_id] = otu
                else:
                    debug("try to add: {} - {}".format(otu_gi_id, otu))
                    debug("but {} was already there".format(self.newadd_gi_otu[otu_gi_id]))
                otu_gi_id = self.data.otu_dict[otu]['^ncbi:gi']
                for item in self.data.otu_dict:
                    if '^ncbi:gi' in self.data.otu_dict[item]:
                        if self.data.otu_dict[item]['^ncbi:gi'] == otu_gi_id and otu != item:
                            exit(-1)
        ###############################
        # debug("self.newseqsgi")
        # debug(self.newseqsgi)
        # debug("self.new_seqs_otu_id")
        # debug(self.new_seqs_otu_id)
        debug("len new seqs dict after remove identical")
        debug(len(self.new_seqs_otu_id))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def find_otudict_gi(self):
        """Used to find seqs that were added twice. Debugging function.
        """
        debug("find_otudict_gi")
        ncbigi_list = []
        for key, val in self.data.otu_dict.items():
            if '^ncbi:gi' in val:
                gi_otu_dict = val["^ncbi:gi"]
                ncbigi_list.append(gi_otu_dict)
        # debug(ncbigi_list)
        return ncbigi_list

    def dump(self, filename=None):
        """writes out class to pickle file

        :param filename: optional filename
        :return: writes out file
        """
        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/scrape_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi)

    def write_query_seqs(self):
        """writes out the query sequence file"""
        debug("write query seq")
        if not self._blast_read:
            self.read_blast()
        self.newseqs_file = "{}.fasta".format(self.date)
        fi = open("{}/{}".format(self.workdir, self.newseqs_file), 'w')
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
            # debug("{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
            # debug(filename)
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        if _VERBOSE:
            sys.stdout.write("aligning query sequences \n")
        # note: sometimes there are still sp in any of the aln/tre
        # hack for the alien taxa thing
        self.remove_alien_aln_tre()
        self.data.write_papara_files()
        os.chdir(self.workdir)  # Clean up dir moving
        try:
            debug("I call papara")
            assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
            # debug(self.newseqs_file)
            subprocess.call(["papara",
                             "-t", "random_resolve.tre",
                             "-s", "aln_ott.phy",
                             # "-j", self.config.num_threads,  # TODO MK: gives error, try to implement for speed up
                             "-q", self.newseqs_file,
                             "-n", papara_runname])  # FIXME directory ugliness
            if _VERBOSE:
                sys.stdout.write("Papara done")
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it installed?\n")
                sys.exit()
            # handle file not found error.
            else:
                # Something else went wrong while trying to run `wget`
                raise
        os.chdir(cwd)
        # debug(os.path.exists(path="{}/papara_alignment.{}".format(self.workdir, papara_runname)))
        assert os.path.exists(path="{}/papara_alignment.{}".format(self.workdir, papara_runname))
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment."
                                                    "{}".format(self.workdir, papara_runname), schema="phylip")
        # debug(self.data.aln.taxon_namespace)
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        # self.data.reconcile()
        # self.data.prune_short()  #
        self._query_seqs_aligned = 1

    def remove_alien_aln_tre(self):
        """Sometimes there were alien entries in self.tre and self.aln.

        This function ensures they are properly removed."""
        for tax_lab in self.data.aln.taxon_namespace:
            if tax_lab not in self.data.tre.taxon_namespace:
                sys.stdout.write("tax not in tre. This is an alien name in the data.")
                # debug("tax not in tre ")
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_lab in self.data.tre.taxon_namespace:
            if tax_lab not in self.data.aln.taxon_namespace:
                sys.stdout.write("tax not in aln. This is an alien name in the data. ")
                self.data.remove_taxa_aln_tre(tax_lab)
        # for tax_aln in self.data.aln.taxon_namespace:
        #     if tax_aln not in self.data.tre.taxon_namespace:
        #         sys.stdout.write("tax not in tre. This is an alien name in the data. ")
        #         self.data.remove_taxa_aln_tre(tax_aln.label)
        # for tax_tre in self.data.tre.taxon_namespace:
        #     if tax_tre not in self.data.aln.taxon_namespace:
        #         sys.stdout.write("tax not in aln. This is an alien name in the data. ")
        #         self.data.remove_taxa_aln_tre(tax_tre.label)
        # added another reconcile step, maybe that helps with the alien taxa problem
        # self.data.reconcile()
        self.data.prune_short()


    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree."""

        # TODO MK: if statement is not yet working, as it contains a path not the number of seqs
        #if len(self.newseqs_file) <= 500:
        if self.backbone is not True:

            if os.path.exists("RAxML_labelledTree.PLACE"):
                os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
            if _VERBOSE:
                sys.stdout.write("placing query sequences \n")
            cwd = (os.getcwd())
            os.chdir(self.workdir)
            try:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
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
                    sys.exit()
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
            os.chdir(cwd)
            self._query_seqs_placed = 1
        else:
            cwd = (os.getcwd())
            os.chdir(self.workdir)
            backbonetre = self.orig_newick
            backbonetre.resolve_polytomies()
            backbonetre.write(path="backbone.tre", schema="newick", unquoted_underscores=True)
            os.chdir(cwd)
            self._query_seqs_placed = 1

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree"""
        cwd = os.getcwd()
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            # debug("{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
            # debug(filename)
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))

        # TODO MK: work on it, first step of not using starting tree was wrong, if that is working un-comment the following stuff
        #if self._query_seqs_placed == 1:  # if too many new sequences are found, do not use a starting tree, is not really faster.
        if self.backbone is not True:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-s", "papara_alignment.extended",
                             "-t", "place_resolve.tre",
                             "-p", "1",
                             "-n", "{}".format(self.date)])
        else:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-s", "papara_alignment.extended",
                             "-r", "backbone.tre",
                             "-p", "1",
                             "-n", "{}".format(self.date)])
        # else:
        #     subprocess.call(["raxmlHPC", "-m", "GTRCAT",
        #                      "-s", "papara_alignment.extended",
        #                      "-p", "1",
        #                      "-n", "{}".format(self.date)])
        os.chdir(cwd)
        self._full_tree_est = 1  # TODO: Currently not used, do we want to use it somewhere?

    def calculate_bootstrap(self):
        """calculate bootstrap and consensus trees
        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        """
        os.chdir(self.workdir)
        # for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
        #     os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        # run bootstrap
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "papara_alignment.extended",
                         # "-t", "place_resolve.tre",
                         "-p", "1", "-b", "1", "-#", "autoMRE",
                         "-n", "{}".format(self.date)])
        # make bipartition tree
        # is the -f b command
        # -z specifies file with multiple trees
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "previous_run/papara_alignment.extended",
                         "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                         "-n", "all{}".format(self.date)])
        # subprocess.call(["raxmlHPC", "-m", "GTRCAT",
        #                  "-s", "previous_run/papara_alignment.extended",
        #                  "-t", "{}/RAxML_bestTree.all{}".format(self.workdir, self.date),
        #                  "-p", "1", "-f", "b", "-z", "RAxML_bootstrap.all{}".format(self.date),
        #                  "-n", "bipart_{}".format(self.date)])
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
        # rapid bootstrapping
        # -f a: is the command to do that
        # subprocess.call(["raxmlHPC", "-m", "GTRCAT",
        #                  "-f", "a",
        #                  "-p", "1", "-x", "1", "-#", "autoMRE", "-s", "aln_ott.phy",
        #                  "-n", "rapidBS_{}".format(self.date)])

    def remove_blacklistitem(self):
        """This removes items from aln, and tree, if the corresponding GI were added to the blacklist.

        Note, that seq that were not added because they were similar to the one being removed here, are lost
        (that should not be a major issue though, as in a new blast_run, they can be added.)
        """
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            acc = self.data.otu_dict[tax.label].get("^ncbi:accession")
            if gi_id in self.blacklist or acc in self.blacklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, gi is part of blacklist"
        # self.data.reconcile()
        self.data.prune_short()

        debug(self.data.tre.as_string(schema='newick'))

    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignment with the expanded ones"""
        debug("generate streamed aln")
        if self.blacklist:
            self.remove_blacklistitem()
        debug(len(self.new_seqs))
        if len(self.new_seqs) > 0:
            # self.remove_identical_seqs()  # speed up
            self.data.write_files()  # should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:  # TODO rename to something more intuitive
                self.write_query_seqs()
                self.align_query_seqs()
                # self.data.reconcile()
                # self.data.prune_short()  # cannot happen here, as aln has new seq but tre not

                self.place_query_seqs()
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
                if self.config.gifilename is not True:
                    os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
                if os.path.exists("{}/last_completed_update".format(self.workdir)):  # TODO: this and the following line are not used.
                    os.rename(self.tmpfi, "{}/last_completed_update".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    # debug(filename.split("/")[-1])
                    if not os.path.exists("{}/previous_run".format(self.workdir)):
                        # debug('{}/previous_run/'.format(self.workdir))
                        os.makedirs('{}/previous_run/'.format(self.workdir))
                    if self.config.gifilename is not True:
                        # if not os.path.exists("{}/previous_run".format(self.workdir)):
                        #     os.makedirs('{}/previous_run/'.format(self.workdir))
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                    else:
                        # debug(filename.split("/")[-1])
                        # debug(self.blast_subdir)
                        # if not os.path.exists("{}/previous_run".format(self.workdir)):
                        #     # debug('{}/previous_run/'.format(self.workdir))
                        #     os.makedirs('{}/previous_run/'.format(self.workdir))
                        # debug( "{}/{}".format(self.workdir, filename.split("/")[-1]))
                        # debug(filename)
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file),
                          "{}/previous_run/newseqs.fasta".format(self.workdir))
                self.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {}  # Wipe for next run
                self.new_seqs_otu_id = {}
                # self.newseqsgi = []  # Never replace it, is used to filter already added gis!!!
                self.repeat = 1
            else:
                if _VERBOSE:
                    sys.stdout.write("No new sequences after filtering.\n")
                self.repeat = 0
        else:
            if _VERBOSE:
                sys.stdout.write("No new sequences found.\n")
            self.repeat = 0
            self.calculate_bootstrap()
        self.reset_markers()
        self.data.dump()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))

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
        # self.sp_d = {}
        # get list of sequences,
        localfiles = os.listdir(path_to_local_seq)
        for index, item in enumerate(localfiles):
            item = str(item)
            if item.startswith(".~"):
                localfiles[index] = None
        localfiles = filter(None, localfiles)
        # debug(localfiles)
        # TODO MK: add assert to ensure that every file is a fasta file in the folder
        for filename in localfiles:
            # debug(file)
            filepath = "{}/{}".format(path_to_local_seq, filename)
            open_file = open(filepath)
            # for line in openFile:
            content = open_file.readlines()
            content = [x.strip() for x in content]
            content = filter(None, content)  # fastest
            count = 0
            gi_id_l = content[::2]
            seq_l = content[1::2]
            # print(gi_id_l, seq_l)
            # in current setup 1 seq per file, but this is written in a way,
            # that a file with multiple seqs can be read in as well
            for i in xrange(0, len(gi_id_l)):
                key = gi_id_l[i].replace(">", "")
                count = count + 1
                seq = seq_l[i]
                # debug(key)
                # debug(seq)
                local_blast.write_blast_files(self.workdir, key, seq, db=True, fn="local_unpubl_seq")
        os.chdir(os.path.join(self.workdir, "blast"))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
        # debug(cmd1)
        os.system(cmd1)

    def make_otu_dict_entry_unpubl(self, key):
        """Adds the local unpublished data to the otu_dict.

        Information are retrieved from the additional json file/self.unpubl_otu_json.
        I make up accession numbers, which might not be necessary

        :param key: unique identifier of the unpublished seq
        :return: generates self.data.gi_dict entry for key
        """
        debug("make_otu_dict_entry_unpubl")
        # debug(self.data.gi_dict.keys())
        # debug(key)
        gi_counter = 1
        if key not in self.data.gi_dict.keys():
            # debug("key is new")
            # numbers starting with 0000 are unpublished data
            self.data.gi_dict[key] = {'accession': "000000{}".format(gi_counter),
                                      'title': "unpublished", 'localID': key[7:]}
            gi_counter += 1
            # self.data.otu_dict[key] = {}
            # self.data.otu_dict[key]['^ncbi:gi'] = key
            # self.data.otu_dict[key]['^ncbi:accession'] = self.data.gi_dict[key]['accession']
            # self.data.otu_dict[key]['^user:TaxonName'] = self.data.gi_dict[key]['localID']
            # self.data.otu_dict[key]['^ncbi:title'] = self.data.gi_dict[key]['title']
            # local_id = self.data.gi_dict[key]['localID']
            # key2 = "otu{}".format(local_id)
            # self.data.otu_dict[key]['^ot:ottTaxonName'] = self.unpubl_otu_json[key2]['^ot:ottTaxonName']
            # self.data.otu_dict[key]['^ncbi:taxon'] = self.unpubl_otu_json[key2]['^ncbi:taxon']
            # self.data.otu_dict[key]['^ot:ottId'] = self.unpubl_otu_json[key2]['^ot:ottId']
            # self.data.otu_dict[key]['^physcraper:status'] = "local seq"
            # self.data.otu_dict[key]['^physcraper:last_blasted'] = "1800/01/01"
            # self.ids.get_rank_info(taxon_name=self.data.otu_dict[key]['^ot:ottTaxonName'])
        else:
            # debug("add new k,v - pairs")
            # debug(self.data.gi_dict[key])
            self.data.gi_dict[key].update([('accession', "000000{}".format(gi_counter)),
                                           ('title', "unpublished"), ('localID', key[7:])])


###############################

class FilterBlast(PhyscraperScrape):
    """Takes the Physcraper Superclass and filters the ncbi blast results to only include a subset of the sequences.

    They can be filtered by number or by rank and number. The feature can be useful for non-population-level studies,
    e.g. analyses which require a single representative per taxon (threshold = 1) or to check the monophyly of taxa
    without having to deal with over-representation of few taxa (e.g. threshold = 4, which allows to get a good overview
    of what is available without having some taxa being represented by high numbers of sequences).
    The second option (downtorank), which is optional, allows to filter according to taxonomic levels, e.g. getting
    a number of representative sequences for a genus or lineage it can also be used to deal with subspecies.

    existing self objects are:
        self.sp_d: dictionary
                key = species name
                value = dictionary:
                    key = otuID
                    value = otu_dict entry
        self.sp_seq_d: dictionary
                key = species name
                value = dictionary (Is overwritten every 'round')
                    key = otuID
                    value = seq.
        self.filtered_seq: dictionary. Is used as the self.new_seqs equivalent from Physcraper, just with fewer seqs.
                            Is overwritten every 'round'
                key = otuID,
                val = seq.
        self.downtorank: optional string defining the level of taxonomic filtering, e.g. "species", "genus"
    """
    # TODO MK: self.sp_d = {} does not need to contain all otu_dict info, key is sufficient

    def __init__(self, data_obj, ids, settings=None):
        super(FilterBlast, self).__init__(data_obj, ids)
        debug("start derived class init")
        # additional things that are needed for the filtering process
        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}
        self.downtorank = None
        # self.localblast = False
        # self.unpubl_otu_json = None
        # self.not_added = []
        # if settings is not None:
        #     self.blacklist = settings.blacklist
        # else:
        #     self.blacklist = []
        # self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,",  "local"]

    def sp_dict(self, downtorank=None):
        """Takes the information from the Physcraper otu_dict and makes a dict with species name as key and
        the corresponding seq information from aln and blast seq, it returns self.sp_d.

        This is generated to make information for the filtering class more easily available. self.sp_d sums up which
        information are available per taxonomic concept and which have not already been removed during either
        the remove_identical_seq steps or during a filtering run of an earlier cycle.

        Note: has test, test_sp_d.py

        :param downtorank: string defining the level of taxonomic filtering, e.g. "species", "genus"
        :return: self.sp_d
        """
        self.downtorank = downtorank
        debug("make sp_dict")
        # debug(self.data.otu_dict)
        self.sp_d = {}
        for key in self.data.otu_dict:
            # debug(key)
            # debug(key['^physcraper:status'])
            if self.data.otu_dict[key]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug(self.downtorank)
                tax_name = self.ids.find_name(sp_dict=self.data.otu_dict[key])
                if tax_name is None:
                    debug("tax_name is None")
                    # debug(self.data.otu_dict[key])
                    gi_id = self.data.otu_dict[key]['^ncbi:gi']
                    if gi_id == None:
                        gi_id = self.data.otu_dict[key]['^ncbi:accession']
                    # debug(gi_id)
                    # debug(type(gi_id))
                    tax_name = self.ids.find_name(gi=gi_id)
                    if tax_name is None:
                        debug("something is going wrong!Check species name")
                        sys.stderr.write("{} has no corresponding tax_name! Check what is wrong!".format(key))
                tax_name = str(tax_name).replace(" ", "_")
                if self.config.blast_loc == 'remote':
                    self.ids.get_rank_info_from_web(taxon_name=tax_name)

                    debug(self.ids.otu_rank.keys())
                    tax_id = self.ids.otu_rank[tax_name]["taxon id"]
                else:
                    tax_id = self.ids.ncbi_parser.get_id_from_name(tax_name)
                if self.downtorank is not None:
                    downtorank_name = None
                    downtorank_id = None
                    tax_name = str(tax_name).replace(" ", "_")
                    if self.config.blast_loc == 'remote':
                        self.ids.get_rank_info_from_web(taxon_name=tax_name)
                        tax_id = self.ids.otu_rank[tax_name]["taxon id"]
                        lineage2ranks = self.ids.otu_rank[str(tax_name).replace(" ", "_")]["rank"]
                        ncbi = NCBITaxa()
                        for key_rank, val in lineage2ranks.iteritems():
                            if val == downtorank:
                                downtorank_id = key_rank
                                value_d = ncbi.get_taxid_translator([downtorank_id])
                                downtorank_name = value_d[int(downtorank_id)]
                    else:
                        tax_id = self.ids.ncbi_parser.get_id_from_name(tax_name)
                        downtorank_id = self.ids.ncbi_parser.get_downtorank_id(tax_id, self.downtorank)
                        downtorank_name = self.ids.ncbi_parser.get_name_from_id(downtorank_id)
                    # following lines should not be necessary, thanks to ncbi_parser
                    # if tax_name not in self.ids.otu_rank.keys():
                    #     self.ids.get_rank_info(taxon_name=tax_name)
                    #     tax_name = str(tax_name).replace(" ", "_")
                    # # if tax_name in self.ids.otu_rank.keys():
                    # lineage2ranks = self.ids.otu_rank[tax_name]["rank"]
                    # # else:
                    # #     self.ids.get_rank_info(taxon_name=tax_name)
                    # #     # debug(self.ids.otu_rank.keys())
                    # #     lineage2ranks = self.ids.otu_rank[str(tax_name).replace(" ", "_")]["rank"]
                    # #     # debug(lineage2ranks)s
                    # ncbi = NCBITaxa()
                    # for key_rank, val in lineage2ranks.iteritems():
                    #     if val == downtorank:
                    #         tax_id = key_rank
                    #         value_d = ncbi.get_taxid_translator([tax_id])
                    #         tax_name = value_d[int(tax_id)]
                    tax_name = downtorank_name
                    tax_id = downtorank_id
                tax_name = tax_name.replace(" ", "_")

                self.ids.spn_to_ncbiid[tax_name] = tax_id
                self.ids.ncbiid_to_spn[tax_id] = tax_name
                # change to tax_id
                if tax_id in self.sp_d:
                    self.sp_d[tax_id].append(self.data.otu_dict[key])
                else:
                    self.sp_d[tax_id] = [self.data.otu_dict[key]]
                # if tax_name in self.sp_d:
                #     self.sp_d[tax_name].append(self.data.otu_dict[key])
                # else:
                #     self.sp_d[tax_name] = [self.data.otu_dict[key]]
        return self.sp_d

    def make_sp_seq_dict(self):
        """Uses the sp_d to make a dict with species names as key1, key2 is gi/sp.name and value is seq

        This is used to select representatives during the filtering step, where it selects how many sequences per
        species to keep in the alignment. It will only contain sp that were not removed in an earlier cycle of the
        program.

        Note: has test, test_sp_seq_d.py

        return: self.sp_seq_d
        """
        debug("make_sp_seq_dict")
        for key in self.sp_d:
            # loop to populate dict. key1 = sp name, key2= gi number, value = seq,
            # number of items in key2 will be filtered according to threshold and already present seq
            # debug(key)
            seq_d = {}
            for otu_id in self.sp_d[key]:
                # following if statement should not be necessary as it is already filtered in the step before.
                # I leave it in for now.
                if '^physcraper:status' in otu_id and otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    # I am using the next if to delimit which seq where already present from an earlier run,
                    # they will get a sp name (str), in order to distinguish them from newly found seq,
                    # which will have the gi (int). This differentiation is needed in the filtering blast step.
                    if otu_id['^physcraper:last_blasted'] != '1800/01/01':
                        tax_name = self.ids.find_name(sp_dict=otu_id)
                        for user_name_aln, seq in self.data.aln.items():
                            otu_dict_label = self.ids.find_name(sp_dict=self.data.otu_dict[user_name_aln.label])
                            if tax_name == otu_dict_label:
                                seq = seq.symbols_as_string().replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[user_name_aln.label] = seq
                    else:
                        if "^ncbi:gi" in otu_id:  # this should not be needed: all new blast seq have gi
                            # debug("gi in otu_id")
                            # gi_num = int(otu_id['^ncbi:gi'])
                            gi_num = otu_id['^ncbi:gi']
                            if gi_num == None:
                                gi_num = otu_id['^ncbi:accession']
                            if gi_num in self.new_seqs.keys():
                                seq = self.new_seqs[gi_num]
                                seq = seq.replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[gi_num] = seq
                    self.sp_seq_d[key] = seq_d
        # debug(self.sp_seq_d)
        return

    def select_seq_by_local_blast(self, seq_d, fn, threshold, count):
        """Selects number of sequences from local_blast to fill up sequences to the threshold. 

        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program 
        how many sequences for the taxon are already available in the aln.

        It will only include species which have a blast score of mean plus/minus sd.
        Uses the information returned by read_local_blast() to select which sequences shall be added in a filtered run.

        Note: has test,test_select_seq_by_local_blast.py

        :param seq_d: is the value of self.sp_d (= another dict)
        :param fn: refers to a filename to find the local blast file produced before,
                    which needs to be read in by read_local_blast()
        :param threshold: threshold
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        # TODO MK: add threshold and downto to self
        debug("select_seq_by_local_blast")
        seq_blast_score = local_blast.read_local_blast(self.workdir, seq_d, fn)
        random_seq_ofsp = {}
        if (threshold - count) <= 0:
            debug("already too many samples of sp in aln, skip adding more.")
        elif len(seq_blast_score.keys()) == (threshold - count):
            random_seq_ofsp = seq_blast_score
        elif len(seq_blast_score.keys()) > (threshold - count):
            random_seq_ofsp = random.sample(seq_blast_score.items(), (threshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        elif len(seq_blast_score.keys()) < (threshold - count):
            random_seq_ofsp = seq_blast_score
        if len(random_seq_ofsp) > 0:
            for key, val in random_seq_ofsp.items():
                self.filtered_seq[key] = val
        # debug(self.filtered_seq)
        return self.filtered_seq

    def select_seq_by_length(self, taxon_id, threshold, count):
        """This is another mode to filter the sequences, if there are more than the threshold.

        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".
        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program how many
        sequences for the taxon are already available in the aln.

        :param taxon_id: key from self.sp_seq_d
        :param threshold: threshold - max number of sequences added per taxon - defined in input
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        debug("select_seq_by_length")
        max_len = max(self.sp_seq_d[taxon_id].values())
        # !!! sometimes the only seq in seq_w_maxlen is the original seq,
        # then this is the one to be added, but it will be removed,
        # later as it is no new seq! thus no new seq for that species is added
        seq_w_maxlen = {}
        for key, val in self.sp_seq_d[taxon_id].iteritems():
            if self.sp_d[taxon_id][key]['^physcraper:status'].split(' ')[0] != ["added", "deleted", "original", "new"]:
                if len(val) == len(max_len):
                    seq_w_maxlen[key] = val
        if (threshold - count) <= 0:
            debug("already to many samples of sp in aln, skip adding more.")
            random_seq_ofsp = None
        elif len(seq_w_maxlen) == (threshold - count):
            random_seq_ofsp = seq_w_maxlen
        elif len(seq_w_maxlen) > (threshold - count):
            random_seq_ofsp = random.sample(seq_w_maxlen.items(), (threshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        else:
            toselect = range(len(seq_w_maxlen), (threshold - count))
            keymax = seq_w_maxlen.keys()
            subdict = {k: v for k, v in self.sp_seq_d[taxon_id].iteritems() if k not in keymax}
            second_len = max(subdict.values())
            seq2len = {}
            for key, val in subdict.iteritems():
                if len(val) == len(second_len):
                    seq2len[key] = val
            random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
            random_seq_ofsp = dict(random_seq_ofsp)
            random_seq_ofsp.update(seq_w_maxlen)
        assert random_seq_ofsp is not None
        if random_seq_ofsp is not None:
            for key in random_seq_ofsp.keys():
                # debug(key)
                self.filtered_seq[key] = random_seq_ofsp[key]

    def add_all(self, key):
        """It adds all seq to filtered_seq dict as the number of sequences present is smaller than the threshold value.

        It is only used, when sequences selection happens via blasting.

        Note: has test, test_add_all.py

        :param key: key of self.sp_d (taxon name)
        :return: self.filtered_seq
        """
        debug('add_all')
        for otu_id in self.sp_d[key]:
            if '^physcraper:status' in otu_id:
                if otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if otu_id['^physcraper:last_blasted'] == '1800/01/01':
                        gi_num = otu_id['^ncbi:gi']
                        if gi_num == None:
                            gi_num = otu_id['^ncbi:accession']
                        # debug(gi_num)
                        seq = self.new_seqs[gi_num]
                        self.filtered_seq[gi_num] = seq
        return self.filtered_seq

    def get_name_for_blastfiles(self, key):
        """Gets the name which is needed to write/read the blast files in 'loop_for_write_blast files'.

        The name needs to be retrieved before the actual loop starts. I use the taxonomic names here,
        as this is the measure of which information goes into which local filter blast database.
        The function is only used within 'loop_for_write_blast_files' to generate the filenames.

        :param key:  key of self.sp_d (taxon name)
        :return: name used for blast file
        """
        nametoreturn = None
        # loop needs to happen before the other one, as we need nametoreturn in second:
        for otu_id in self.sp_d[key]:
            # this if should not be necessary
            tax_name = self.ids.find_name(sp_dict=otu_id)
            if '^physcraper:status' in otu_id and otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if otu_id['^physcraper:last_blasted'] != '1800/01/01':
                # if tax_name is None:
                #     # debug(key)
                #     tax_name = self.ids.get_rank_info(taxon_name=key)
                #     # otu_id['^ot:ottTaxonName'] = tax_name
                # tax_name = tax_name.replace(" ", "_")
                for tax_name_aln, seq in self.data.aln.items():
                    otu_dict_name = self.ids.find_name(sp_dict=self.data.otu_dict[tax_name_aln.label])
                    if tax_name == otu_dict_name:
                        nametoreturn = tax_name_aln.label
            # the next lines where added because the test was breaking,
            # needs thorough testing if it not breaks something else now.
            assert tax_name is not None  # assert instead of if
            if nametoreturn is None and tax_name is not None:
                nametoreturn = tax_name.replace(" ", "_")
            # debug(nametoreturn)
            # else:
            #     debug("do something?")
        # return nametoreturn
        return key

    def loop_for_write_blast_files(self, key):
        """This loop is needed to be able to write the local blast files for the filtering step correctly.

        Function returns a filename for the filter blast, which were generated with 'get_name_for_blastfiles()'.

        Note: has test,test_loop_for_blast.py

        :param key: key of self.sp_d (taxon name)
        :return: name of the blast file
        """
        debug("length of sp_d key")
        # debug(len(self.sp_d[key]))
        # debug(key)
        nametoreturn = self.get_name_for_blastfiles(key)
        for otu_id in self.sp_d[key]:
            # debug("in writing file for-loop")
            # this if should not be necessary, I leave it in for now
            if '^physcraper:status' in otu_id and otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug("generate files used for blast")
                if otu_id['^physcraper:last_blasted'] != '1800/01/01':  # old seq
                    tax_name = self.ids.find_name(sp_dict=otu_id)
                    for tax_name_aln, seq in self.data.aln.items():
                        otu_dict_name = self.ids.find_name(sp_dict=self.data.otu_dict[tax_name_aln.label])
                        if tax_name == otu_dict_name:
                            filename = nametoreturn
                            # filename = tax_name_aln.label
                            if self.downtorank is not None:
                                filename = key
                            # print(filename, seq)
                            local_blast.write_blast_files(self.workdir, filename, seq)
                else:
                    # debug("make gilist as local blast database")
                    if "^ncbi:gi" in otu_id:
                        gi_num = int(otu_id['^ncbi:gi'])
                        if gi_num == None:
                            gi_num = otu_id['^ncbi:accession']
                        # debug(gi_num)
                        # debug("new")
                        file_present = False
                        # debug(gi_num)
                        if gi_num in self.new_seqs.keys():
                            file_present = True
                        if file_present:  # short for if file_present == True
                            if '^physcraper:status' in otu_id:
                                if otu_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                                    filename = gi_num
                                    # debug("write seq to db")
                                    # debug(nametoreturn)
                                    seq = self.sp_seq_d[key][gi_num]
                                    if self.downtorank is not None:
                                        filename = key
                                        nametoreturn = key
                                    # debug(filename)
                                    local_blast.write_blast_files(self.workdir, filename, seq, db=True, fn=nametoreturn)
                                    # blastfile_taxon_names[gi_num] = gi_num
                    namegi = key
        if self.downtorank is not None:
            nametoreturn = key
        if nametoreturn is None:
            nametoreturn = namegi
        assert nametoreturn is not None
        return nametoreturn

    def count_num_seq(self, tax_id):
        """Counts how many sequences there are for a tax_name, excluding sequences that have not been added
        during earlier cycles.

        Function is only used in how_many_sp_to_keep().

        :param tax_id: key from self.sp_seq_d
        :return: dict which contains information of how many seq are already present in aln, how many new seq have been
                found and if the taxon is a new taxon or if seq are already present
        """
        seq_present = 0
        if tax_id in self.sp_seq_d.keys():
            for sp_keys in self.sp_seq_d[tax_id].keys():
                if isinstance(sp_keys, str):
                    seq_present += 1
                if isinstance(sp_keys, unicode):
                    seq_present += 1
        # this determines if a taxonomic name / otu is already present in the aln and how many new seqs were found
        new_taxon = True
        query_count = 0
        for item in self.sp_d[tax_id]:
            if '^physcraper:status' in item and item['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                if item['^physcraper:last_blasted'] != '1800/01/01':
                    new_taxon = False
                if item['^physcraper:status'] == "query":
                    query_count += 1
        count_dict = {'seq_present': seq_present, 'query_count': query_count, 'new_taxon': new_taxon}
        return count_dict

    def how_many_sp_to_keep(self, threshold, selectby):
        """Uses the sp_seq_d and places the number of sequences according to threshold into the self.filterdseq_dict.

        This is essentially the key function of the Filter-class, it wraps up everything.

        :param threshold: threshold value, defined in input
        :param selectby: mode of sequence selection, defined in input
        :return: nothing specific, it is the function, that completes the self.filtered_seq, which contains the filtered
                sequences that shall be added.
        """
        debug("how_many_sp_to_keep")
        debug("length of sp_d")
        debug(len(self.sp_d))
        for tax_id in self.sp_d:
            debug(tax_id)
            count_dict = self.count_num_seq(tax_id)
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            if len(self.sp_d[tax_id]) <= threshold:  # add all stuff to self.filtered_seq[gi_n]
                self.add_all(tax_id)
            elif len(self.sp_d[tax_id]) > threshold:  # filter number of sequences
                debug("filter number of sequences")
                # debug(self.sp_seq_d[tax_id].keys())
                if tax_id in self.sp_seq_d.keys():
                    if selectby == "length":
                        # debug("{}, {}, {}".format(self.sp_seq_d[tax_id], threshold, seq_present))
                        self.select_seq_by_length(tax_id, threshold, seq_present)
#                        self.select_seq_by_length(self.sp_seq_d[tax_id], threshold, seq_present)
                    elif selectby == "blast":
                        # if seq_present >= 1 and seq_present < threshold and count_dict["new_taxon"] == False and query_count != 0:
                        if 1 <= seq_present < threshold and count_dict["new_taxon"] is False and query_count != 0:
                            # debug("seq_present>0")
                            # species is not new in alignment, make blast with existing seq
                            if query_count + seq_present > threshold:
                                taxonfn = self.loop_for_write_blast_files(tax_id)
                                # # next loop does not seem to be used
                                # for element in self.sp_d[tax_id]:
                                #     if '^ot:ottTaxonName' in element:
                                #         blast_seq = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                #         blast_db = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                if self.downtorank is not None:
                                    taxonfn = tax_id
                                local_blast.run_local_blast(self.workdir, taxonfn, taxonfn)
                                self.select_seq_by_local_blast(self.sp_seq_d[tax_id], taxonfn, threshold, seq_present)
                            elif query_count + seq_present <= threshold:
                                self.add_all(tax_id)
                        # species is completely new in alignment
                        elif seq_present == 0 and count_dict["new_taxon"] is True and query_count >= 1:
                            # debug("completely new taxon to blast")
                            # species is completely new in alignment, make blast with random species
                            # debug(count_dict)
                            # debug(tax_id)
                            # debug(self.sp_seq_d)
                            # this causes to add some taxa twice to aln and phy!!! never use add_otu twice!
                            # for item in self.sp_d[tax_id]:
                            #     if '^ncbi:gi' in item:
                            #         # if self.config.blast_loc == 'local':
                            #         #     localblast = True
                            #         # else:
                            #         #     localblast = False
                            #         self.data.add_otu(item['^ncbi:gi'], self.ids)
                            blast_seq = self.sp_seq_d[tax_id].keys()[0]
                            if self.downtorank is not None:
                                str_db = tax_id
                            else:
                                if type(blast_seq) == int:
                                    str_db = str(tax_id)
                                else:
                                    str_db = str(blast_seq)
                            # write files for local blast first:
                            seq = self.sp_seq_d[tax_id][blast_seq]
                            local_blast.write_blast_files(self.workdir, str_db, seq)  # blast qguy
                            # debug("blast db new")
                            blast_db = self.sp_seq_d[tax_id].keys()[1:]
                            # debug(blast_db)
                            for blast_key in blast_db:
                                seq = self.sp_seq_d[tax_id][blast_key]
                                local_blast.write_blast_files(self.workdir, blast_key, seq, db=True, fn=str_db)
                            # make local blast of sequences
                            local_blast.run_local_blast(self.workdir, str_db, str_db)
                            if len(self.sp_seq_d[tax_id]) + seq_present >= threshold:
                                self.select_seq_by_local_blast(self.sp_seq_d[tax_id], str_db, threshold, seq_present)
                            elif len(self.sp_seq_d[tax_id]) + seq_present < threshold:
                                self.add_all(tax_id)
                else:
                    debug("taxon not in sp_seq_dict")
        return

    def replace_new_seq(self):
        """Function to replace self.new_seqs and self.new_seqs_otu_id with the subset of filtered sequences.

        This is the final step in the FilterBlast class, from here it goes back to PhyScraper.

        :return: subsets of self.new_seqs and self.new_seqs_otu_id
        """
        debug("replace new seq")
        # debug(self.filtered_seq)
        keylist = self.filtered_seq.keys()
        # debug(keylist)
        if not self.unpublished:
            keylist = [x for x in keylist if type(x) == int]
        # debug(self.new_seqs.keys())
        seq_not_added = self.new_seqs.keys()
        seq_not_added = [x for x in seq_not_added if type(x) == int]
        reduced_new_seqs_dic = {}
        for gi_id in seq_not_added:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_id:
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'not added, there are enough seq per sp in tre'
        for gi_id in keylist:
            for key in self.data.otu_dict.keys():
                # debug(self.data.otu_dict[key])
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_id:
                        reduced_new_seqs_dic[key] = self.filtered_seq[gi_id]
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added, as representative of taxon'
                        # self.data.otu_dict[otu_id]['^ncbi:gi'] = gi_id
        # might not be necessary, I was missing some gi's when I was replacing the original one.
        # i leave the code here for now
        # reduced_gi_dict = {k: self.data.gi_dict[k] for k in keylist}
        # debug(reduced_gi_dict)
        # self.data.gi_dict.clear()
        # self.data.gi_dict = reduced_gi_dict # data.gi_dict seems to only have newly blasted stuff
        reduced_new_seqs = {k: self.filtered_seq[k] for k in keylist}
        # debug(reduced_new_seqs_dic)
        with open(self.logfile, "a") as log:
            log.write("{} sequences added after filtering, of {} before filtering\n".format(len(reduced_new_seqs_dic),
                                                                                            len(self.new_seqs_otu_id)))
        # debug("self.new_seqs")
        # debug(self.new_seqs)
        # debug(len(self.new_seqs))
        # debug("reduced_new_seqs")
        # debug(reduced_new_seqs)
        # debug(len(reduced_new_seqs))
        # debug("self.new_seqs")
        # debug(self.new_seqs)
        # debug(len(self.new_seqs))
        self.new_seqs = deepcopy(reduced_new_seqs)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic)
        # set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return

    def write_otu_info(self, downtorank=None):
        """Writes different output tables to file: Makes reading important information less code heavy.

        1. table with taxon names and sampling.
        2. a file with all relevant GenBank info to file (otu_dict).

        It uses the self.sp_d to get sampling information, that's why the downtorank is required.

        :param downtorank: hierarchical filter
        :return: writes output to file
        """
        debug("write out infos")
        if len(self.sp_d) == 0:
            sp_d = self.sp_dict(downtorank)
        else:
            sp_d = self.sp_d
        sp_info = {}
        for k in sp_d:
            sp_info[k] = len(sp_d[k])
        with open('taxon_sampling.csv', 'w') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in sp_info.items():
                writer.writerow([key, value])
        otu_dict_keys = ['^ot:ottTaxonName', '^ncbi:gi', '^ncbi:accession', '^physcraper:last_blasted',
                         '^physcraper:status', '^ot:ottId', '^ncbi:taxon', '^ncbi:title']
        with open('otu_seq_info.csv', 'w') as output:
            writer = csv.writer(output)
            for otu in self.data.otu_dict.keys():
                rowinfo = [otu]
                for item in otu_dict_keys:
                    if item in self.data.otu_dict[otu].keys():
                        tofile = str(self.data.otu_dict[otu][item]).replace("_", " ")
                        rowinfo.append(tofile)
                    else:
                        rowinfo.append("-")
                writer.writerow(rowinfo)

    def print_sp_d_as_is(self):
        if self.sp_d is not None:
            with open('sp_d_info.csv', 'w') as output:
                writer = csv.writer(output)
                for group in self.sp_d.keys():
                    rowinfo = [group]
                    for otu in self.sp_d[group]:
                        for item in otu.keys():
                            tofile = str(otu[item]).replace("_", " ")
                            rowinfo.append(tofile)
                    writer.writerow(rowinfo)

    def print_sp_d_recalc(self, downtorank):
        if self.sp_d is not None:
            sp_d = self.sp_dict(downtorank)
            with open('sp_d_info_recalc.csv', 'w') as output:
                writer = csv.writer(output)
                for group in sp_d.keys():
                    rowinfo = list(str(group))
                    debug(rowinfo)
                    debug(sp_d[group])
                    for otu in sp_d[group]:
                        debug(otu)
                        # debug(sp_d[group][otu])
                        # debug(sp_d[group][otu].keys())
                        for item in otu.keys():
                            tofile = str(otu[item]).replace("_", " ")
                            debug(tofile)
                            # tofile = str(sp_d[group][otu][item]).replace("_", " ")
                            rowinfo.append(otu)
                            rowinfo.append(tofile)
                    writer.writerow(rowinfo)


class Settings(object):
    """A class to store all settings for PhyScraper.
    """

    def __init__(self, seqaln, mattype, trfn, schema_trf, workdir, threshold=None,
                 selectby=None, downtorank=None, spInfoDict=None, add_unpubl_seq=None,
                 id_to_spn_addseq_json=None, configfi=None, blacklist=None, shared_blast_folder=None,
                 delay=None, trim=None):

        """Initialize the settings."""
        self.seqaln = seqaln
        self.mattype = mattype
        self.trfn = trfn
        self.schema_trf = schema_trf
        self.workdir = workdir
        self.threshold = threshold
        self.selectby = selectby
        self.downtorank = downtorank
        self.spInfoDict = spInfoDict
        self.add_unpubl_seq = add_unpubl_seq
        self.id_to_spn_addseq_json = id_to_spn_addseq_json
        self.configfi = configfi
        self.blacklist = blacklist
        self.shared_blast_folder = shared_blast_folder
        self.delay = delay
        self.trim = trim


####################


def get_ncbi_tax_id(handle):
    """Get the taxon ID from ncbi.

    :param handle: NCBI read.handle
    :return: ncbi_id
    """
    ncbi_id = None
    gb_list = handle[0]['GBSeq_feature-table'][0]['GBFeature_quals']
    for item in gb_list:
        if item[u'GBQualifier_name'] == 'db_xref':
            # debug(item[u'GBQualifier_value'])
            if item[u'GBQualifier_value'][:5] == 'taxon':
                ncbi_id = int(item[u'GBQualifier_value'][6:])
                break
            else:
                continue
    return ncbi_id


def get_ncbi_tax_name(handle):
    """Get the sp name from ncbi.

    :param handle: NCBI read.handle
    :return: ncbi_spn
    """
    ncbi_sp = None
    gb_list = handle[0]['GBSeq_feature-table'][0]['GBFeature_quals']
    for item in gb_list:
        if item[u'GBQualifier_name'] == 'organism':
            ncbi_sp = str(item[u'GBQualifier_value'])
            ncbi_sp = ncbi_sp.replace(" ", "_")
    return ncbi_sp



# from guppy import hpy

# h = hpy()

# print h.heap() 
