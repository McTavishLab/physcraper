#!/usr/bin/env python
"""Physcraper module"""
import sys
import re
import os
import csv
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import urllib2
from copy import deepcopy
from ete2 import NCBITaxa
import physcraper.AWSWWW as AWSWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from dendropy import Tree, \
    DnaCharacterMatrix, \
    DataSet, \
    datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI, APIWrapper
from peyotl.sugar import tree_of_life, taxomachine
from peyotl.nexson_syntax import extract_tree, \
    get_subtree_otus, \
    extract_otu_nexson, \
    PhyloSchema
# extension functions
import concat  # is the local concat class
import ncbi_data_parser  # is the ncbi data parser class and associated functions
import local_blast

_DEBUG = 0
_DEBUG_MK = 0

_VERBOSE = 0


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


# which python physcraper file do I use?
debug("Current --init-- version number: 10-15-2018.0")
debug(os.path.realpath(__file__))


def get_raw_input():
    """Asks for yes or no user input.

    :return: user input
    """
    debug("get raw input")
    is_valid = 0
    x = None
    while not is_valid:
        try:
            x = raw_input("Please write either 'yes' or 'no': ")
            if x == "yes" or x == "no":
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError, e:
            print ("'%s' is not a valid answer." % e.args[0].split(": ")[1])
    return x


class ConfigObj(object):
    """
    #### To build the class the following is needed:
      * **configfi**: a configuration file in a specific format,
       
            e.g. to read in self.e_value_thresh.
                The file needs to have a heading of the format: [blast] and then somewhere below that heading
                a string e_value_thresh = value

      * **interactive**: defaults to True, is used to interactively update the local blast databases

    #### During the initializing process the following self objects are generated:
      * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out:
                                https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
      * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
      * **self.seq_len_perc**: value from 0 to 1. Defines how much shorter new seq can be compared to input
      * **self.get_ncbi_taxonomy**: Path to sh file doing something...
      * **self.ncbi_dmp**: path to file that has gi numbers and the corresponding ncbi tax id's
      * **self.phylesystem_loc**: defines which phylesystem for OpenTree datastore is used. 
                                The default is api, but can run on local version too. 
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
      * **optional self.objects**:
          * if blastloc == local:
              * self.blastdb: this defines the path to the local blast database
              * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
              * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
    """

    def __init__(self, configfi, interactive=True):
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
        self.phylesystem_loc = config['phylesystem']['location']
        assert (self.phylesystem_loc in ['local', 'api'])  # default is api, but can run on local version of OpenTree datastore
        self.ott_ncbi = config['taxonomy']['ott_ncbi']  # TODO: how do we update the file?
        assert os.path.isfile(self.ott_ncbi)
        # rewrites relative path to absolute path so that it behaves when changing dirs
        self.id_pickle = os.path.abspath(config['taxonomy']['id_pickle'])
        self.email = config['blast']['Entrez.email']
        assert '@' in self.email
        self.blast_loc = config['blast']['location']
        self.num_threads = config['blast'].get('num_threads')
        assert self.blast_loc in ['local', 'remote']
        if self.blast_loc == 'local':
            self.blastdb = config['blast']['localblastdb']
            self.url_base = None
            self.ncbi_parser_nodes_fn = config['ncbi_parser']["nodes_fn"]
            self.ncbi_parser_names_fn = config['ncbi_parser']["names_fn"]
        if self.blast_loc == 'remote':
            self.url_base = config['blast'].get('url_base')
        self.gb_id_filename = config['blast'].get('gb_id_filename', False)
        if self.gb_id_filename is not False:
            if self.gb_id_filename == "True" or self.gb_id_filename == "true":
                self.gb_id_filename = True
            else:
                self.gb_id_filename = False
        if interactive is True:
            self._download_ncbi_parser()
            self._download_localblastdb()
        self.unmapped = config['blast']['unmapped']
        assert self.unmapped in ['remove', 'keep']

        debug("shared blast folder?")
        debug(self.gb_id_filename)
        debug("check db file status?")
        debug(interactive)
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
                    os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'" 
                              "{}/".format(self.blastdb))
                    os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'" 
                              "{}/".format(self.blastdb))
                    cwd = os.getcwd()
                    os.chdir(self.blastdb)
                    os.system("update_blastdb nt")
                    os.system("cat *.tar.gz | tar -xvzf - -i")
                    os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                    os.chdir(cwd)
                elif x == "no":
                    print("You did not agree to download data from ncbi. Program will default to blast web-queries.")
                    print("This is slow and crashes regularly!")
                    self.blast_loc = 'remote'
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime("{}/nt.01.nhr".format(self.blastdb))
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days
                if time_passed >= 90:
                    print("Your databases might not be uptodate anymore. You downloaded them {} days ago. "
                          "Do you want to update the blast databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms".format(time_passed))
                    x = get_raw_input()
                    if x == "yes":
                        cwd = os.getcwd()
                        os.chdir(self.blastdb)
                        os.system('update_blastdb nt')
                        os.system('cat *.tar.gz | tar -xvzf - -i')
                        os.system("update_blastdb taxdb")
                        os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                        os.chdir(cwd)
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type 'yes' or 'no'!")

    def _download_ncbi_parser(self):
        """Check if files are present and if they are up to date.
        If not files will be downloaded. 
        """
        if self.blast_loc == 'local':
            if not os.path.isfile(self.ncbi_parser_nodes_fn):
                print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms")
                x = get_raw_input()
                if x == "yes":
                    os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                    os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                    os.system("mv nodes.dmp ./tests/data/")
                    os.system("mv names.dmp ./tests/data/")
                elif x == "no":
                    print("You did not agree to download data from ncbi. Program will default to blast web-queries.")
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
                if time_passed >= 90:
                    print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms")
                    x = get_raw_input()
                    if x == "yes":
                        os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                        os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                        os.system("mv nodes.dmp ./tests/data/")
                        os.system("mv names.dmp ./tests/data/")
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
    except urllib2.HTTPError, err:
        sys.stderr.write(err)
    # except:  # TODO: seems to be an http error. Did not fgure out how to handle them (requests.exceptions.HTTPError)
        sys.stderr.write("couldn't find study id {} in phylesystem location {}\n".format(study_id, phylesystem_loc))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treebase record associated with study ")
        sys.exit()
    else:
        tb_id = treebase_url.split(':S')[1]
        url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexus".format(tb_id)
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

    Study and tree ID's can be obtained by using scripts/find_trees.py LINEAGE_NAME

    Input can be either a study ID and tree ID from OpenTree
    # TODO: According to code it cannot be either, but must be both

    Alignment need to be a Dendropy DNA character matrix!

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
    # this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    otus = get_subtree_otus(nexson, tree_id=tree_id)
    otu_dict = {}
    orig_lab_to_otu = {}
    treed_taxa = {}
    for otu_id in otus:
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
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)
    # newick should be bare, but alignment should be DNACharacterMatrix


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
    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        otu_dict = json.load(open(otu_json, "r"))
        ott_ids = [otu_dict[otu].get(u'^ot:ottId', ) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)


def standardize_label(item):
    """Make sure that the tipnames are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

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

    :param spp_name: species name
    :return:
    """
    debug(spp_name)
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
    with open(id_to_spn, mode='r') as infile:
        for lin in infile:
            ottid, ottname, ncbiid = None, None, None
            tipname, species = lin.strip().split(',')
            clean_lab = standardize_label(tipname)
            assert clean_lab not in sp_info_dict
            otu_id = "otu{}".format(clean_lab)
            spn = species.replace("_", " ")
            info = get_ott_taxon_info(spn)
            if info:
                ottid, ottname, ncbiid = info
            if not info:
                if _DEBUG:
                    sys.stdout.write("match to taxon {} not found in open tree taxonomy. Trying NCBI next.\n".format(spn))
                ncbi = NCBITaxa()
                name2taxid = ncbi.get_name_translator([spn])
                if len(name2taxid.items()) >= 1:
                    if _DEBUG:
                        sys.stdout.write("found taxon {} in ncbi".format(spn))
                    ncbiid = name2taxid.items()[0][1][0]
                    ottid = id_dict.ncbi_to_ott[ncbiid]
                    ottname = id_dict.ott_to_name[ottid]
                else:
                    sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                     "Proceeding without taxon info\n".format(spn))
            sp_info_dict[otu_id] = {'^ncbi:taxon': ncbiid, '^ot:ottTaxonName': ottname, '^ot:ottId': ottid,
                                    '^ot:originalLabel': tipname, '^user:TaxonName': species,
                                    '^physcraper:status': 'original', '^physcraper:last_blasted': "1900/01/01"}
    return sp_info_dict


class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
        Hypothetically, all the keys in the  otu_dict should be clean.

        ####To build the class the following is needed:

          * **newick**: dendropy.tre.as_string(schema=schema_trf) object
          * **otu_dict**: json file including the otu_dict information generated earlier
          * **alignment**: dendropy DNACharacterMatrix object
          * **ingroup_mrca**: OToL identifier of the group of interest, either subclade as defined by user or of
                            all tiplabels in the phylogeny
          * **workdir**: the path to the corresponding working directory
          * **schema**: optional argument to define tre file schema, if different from "newick"

        ####During the initializing process the following self objects are generated:

          * **self.aln**: contains the alignment and which will be updated during the run
          * **self.tre**: contains the phylogeny, which will be updated during the run
          * **self.otu_dict**: dictionary with taxon information and physcraper relevant stuff
               * key: a unique identifier (otu plus either "tiplabel of phylogeny" or for newly found sequences
                    PS_number.
               * value: dictionary with the following key:values:
                    * '^ncbi:gi': GenBank identifier - deprecated by Genbank - only older sequences will have it
                    * '^ncbi:accession': Genbanks accession number
                    * '^ncbi:title': title of Genbank sequence submission
                    * '^ncbi:taxon': ncbi taxon identifier
                    * '^ot:ottId': OToL taxon identifier
                    * '^physcraper:status': contains information if it was 'original', 'queried', 'removed',
                                        'added during filtering process'
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

        ####Following functions are called during the init-process:

          * **self._reconcile_names()**: 
                removes taxa, that are not found in both, the phylogeny and the aln and changes their names????

        ####The physcraper class is then updating: 
          * self.aln, self.tre and self.otu_dict, self.ps_otu, self.gi_dict
 """

    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None, taxon_namespace=None):
        # TODO add assertions that inputs are correct type!!!
        debug("build ATT class")
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
        assert (self.tre.taxon_namespace is self.aln.taxon_namespace)
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
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
        self.gb_dict = {}  # has all info about new blast seq TODO: Cannot be deleted. Is used frequently! TODODELTE (should maybe go anyhow due to gi switch?): Should this not be part of physcraper class instead? it has all blast information. Blast is not part of this class.
        self._reconciled = False  # TODO: for what do we want to use it? .... it was checking to see if name reconcilation has ahappened yet. Should get flipped to true when done. MK: Yes, but we never do anything with the information
        self.unpubl_otu_json = None

    def _reconcile_names(self):
        """Taxa that are only found in the tree, or only in the alignment are deleted.

        This checks that the tree "original labels" from phylsystem
        align with those found in the alignment. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!
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
            # potentially slow at large number of taxa and large numbers to be pruned
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
                    original = self.otu_dict[otu].get('^ot:originalLabel')
                    if original == tax.label or original == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))
    # TODO - make sure all taxon labels are unique OTU ids.

    def prune_short(self, min_seqlen_perc=0.75):
        """Prunes sequences from alignment if they are shorter than 75%, or if tip is only present in tre.

        Sometimes in the de-concatenating of the original alignment taxa with no sequence are generated
        or in general if certain sequences are really short. This removes those from both the tre and the alignment.

        has test: test_prune_short.py

        :param min_seqlen_perc: minimum length of seq
        :return: prunes aln and tre
        """
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
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short "
                     "step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                self.remove_taxa_aln_tre(tax.label)
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in prune short"
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        assert treed_taxa.issubset(aln_ids)
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        self.trim()
        self._reconciled = 1

    def trim(self, taxon_missingness=0.75):
        """ It removes bases at the start and end of alignments, if they are represented by less than 75%
        of the sequences in the alignment.
        Ensures, that not whole chromosomes get dragged in. It's cutting the ends of long sequences.

        Used in prune_short()
        has test: test_trim.py

        :param taxon_missingness: defines how many sequences need to have a base at the start/end of an alignment
        """
        # debug('in trim')
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
            aln_ids.add(taxon.label)
        assert aln_ids.issubset(self.otu_dict.keys())
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        for leaf in self.tre.leaf_nodes():
            if leaf.taxon not in aln_ids:
                self.tre.prune_taxa([leaf])
                self.tre.prune_taxa_with_labels([leaf.taxon])
                self.tre.prune_taxa_with_labels([leaf])
                treed_taxa.remove(leaf.taxon)
        assert treed_taxa.issubset(aln_ids)
        if _VERBOSE:
            sys.stdout.write("trimmed alignment ends to < {} missing taxa, "
                             "start {}, stop {}\n".format(taxon_missingness, start, stop))
        return

    def add_otu(self, gi, ids_obj):
        """generates an otu_id for new sequences and adds them into the otu_dict.
        Needs to be passed an IdDict to do the mapping"""
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        ncbi_id = ids_obj.map_gi_ncbi(gi)  # check that try/except is working here
        # TODO do we need rank info here?

        # seems like without the try and except we are missing tons of information
        try:
            # ncbi_id = int(ids_obj.map_gi_ncbi(gi))
            # try:
            ott = int(ids_obj.ncbi_to_ott[ncbi_id])
            # except:
            if ott is None:
                ott = "OTT_{}".format(self.ps_otu)
                self.ps_otu += 1
            spn = str(ids_obj.ott_to_name[ott]).replace(" ", "_")
        except:
            spn = ids_obj.get_rank_info(gi_id=gi)

            ncbi_id = ids_obj.otu_rank[spn]["taxon id"]
            try:
                ott = int(ids_obj.ncbi_to_ott[ncbi_id])
            except:
                ott = "OTT_{}".format(self.ps_otu)
                self.ps_otu += 1
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi]['title']
        self.otu_dict[otu_id]['^ncbi:taxon'] = ncbi_id
        self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ncbi_id)
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = ids_obj.ott_to_name.get(self.otu_dict[otu_id]['^ot:ottId'])
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"  # 1800 = never blasted; 1900 = blasted 1x, not added; this century = blasted and added
        if _DEBUG >= 2:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi, otu_id))
        return otu_id

    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """Papara is finicky about trees and needs phylip, this writes out needed files for papara
        (except query sequences)"""
        # CAN I even evaulte things in the function definitions?
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
        # First write rich annotation json file with everything needed for later?
        debug("write_files")
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
        :param treepath: optional: full filenname (including path) for phylogeny
        :param alnpath:  optional: full filenname (including path) for alignment
        :param norepeats: optional: if there shall be no duplicate names in the labelled output files
        :param add_gb_id: optional, to supplement tiplabel with corresponding GenBank sequence identifier
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
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label is None:
                if self.otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_{}".format(self.otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"), self.otu_dict[taxon.label].get('^ot:ottTaxonName', "unk"))
            new_label = str(new_label).replace(' ', '_')
            if gi_id:
                sp_counter = 2
                if new_label in new_names and norepeats:
                    gi_id = self.otu_dict[taxon.label].get('^ncbi:gi')
                    if gi_id is None:
                        gi_id = sp_counter
                        sp_counter += 1
                    new_label = "_".join([new_label, str(gi_id)])
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
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found
            # the first prune does not remove it sometimes...
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, updated otu_dict but was never in tre or aln!"

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
    if None in ott_ids:
        ott_ids.remove(None)
    synth_tree_ott_ids = []
    ott_ids_not_in_synth = []
    for ott in ott_ids:
        try:
            tree_of_life.mrca(ott_ids=[ott], wrap_response=False)
            synth_tree_ott_ids.append(ott)
        except:
            ott_ids_not_in_synth.append(ott)
    if len(synth_tree_ott_ids) == 0:
        sys.stderr.write('No sampled taxa were found in the current sysnthetic tree. Please find and input and approppriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit()
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
        print(mrca_node.keys())
        sys.stderr.write('(v3) MRCA of sampled taxa not found. Please find and input and approppriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit()
    return tax_id


def get_ott_ids_from_otu_dict(otu_dict):  # TODO put into data obj?
    """Get the ott ids from an otu dict object"""
    ott_ids = []
    for otu in otu_dict:
        try:
            ott_ids.append(otu['^ot:ottId'])
        except KeyError:
            pass


#####################################

class IdDicts(object):
    """Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

        ####To build the class the following is needed:
            
          * **config_obj**: Object of class config (see above)
          * **workdir**: the path to the assigned working directory
          * **mrca**: mrca as defined by input, can be a class

        ####During the initializing process the following self objects are generated:
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

          * **self.mrca_ott**: user defined list of mrca OTT-ID's
          * **self.mrca_ncbi**: set, which is fed by self.get_ncbi_mrca()

          * **Optional**:
              * depending on blasting method:
               * self.ncbi_parser: for local blast, initializes the ncbi_parser class, that contains information about rank and identifiers
               * self.otu_rank: for remote blast to store the rank information
    """

    # TODO - could - should be shared acrosss runs?! .... nooo.
    def __init__(self, config_obj, workdir):
        """Generates a series of name disambiguation dicts"""
        self.workdir = workdir
        self.config = config_obj
        assert self.config.email
        self.ott_to_ncbi = {}
        self.ncbi_to_ott = {}
        self.ott_to_name = {}
        self.gi_ncbi_dict = {}
        self.otu_rank = {}
        fi = open(config_obj.ott_ncbi)  # TODO need to keep updated
        for lin in fi:  # TODO This is insanely memory inefficient
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()
            assert len(self.ott_to_ncbi) > 0
            assert len(self.ncbi_to_ott) > 0
            assert len(self.ott_to_name) > 0
        fi.close()
        if os.path.isfile("{}/id_map.txt".format(workdir)):  # todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]

    def get_rank_info(self, gi_id=False, taxon_name=False):
        """Collects rank and lineage information from ncbi,
        used to delimit the sequences from blast,
        when you have a local blast database or a Filter Blast run
        """
        # debug("get_rank_info")
        Entrez.email = self.config.email
        if gi_id:
            # debug("gi_id to tax_name")
            tries = 5
            for i in range(tries):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
                except:
                    if i < tries - 1:  # i is zero indexed
                        continue
                    else:
                        raise
                break
            read_handle = Entrez.read(handle)[0]
            tax_name = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value']
        else:
            tax_name = str(taxon_name).replace("_", " ")
        if tax_name not in self.otu_rank.keys():
            # debug("tax_name to rank")
            ncbi = NCBITaxa()
            try:
                tax_id = int(Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0])
            except:
                # debug("except")
                tax_info = ncbi.get_name_translator([tax_name])
                if tax_info == {}:
                    print("Taxon name does not match any species name in ncbi. Check that name is written correctly!")
                tax_id = int(tax_info.items()[0][1][0])
            ncbi = NCBITaxa()
            lineage = ncbi.get_lineage(tax_id)
            lineage2ranks = ncbi.get_rank(lineage)
            tax_name = str(tax_name).replace(" ", "_")
            assert type(tax_id) == int
            self.otu_rank[tax_name] = {"taxon id": tax_id, "lineage": lineage, "rank": lineage2ranks}
        return tax_name

    def map_gi_ncbi(self, gi):
        """get the ncbi taxon id's for a gi input"""
        if _DEBUG == 2:
            sys.stderr.write("mapping gi {}\n".format(gi))
        tax_id = None
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
            tax_name = self.get_rank_info(gi_id=gi)
            tax_id = self.otu_rank[tax_name]["taxon id"]
            self.gi_ncbi_dict[gi] = tax_id
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
    """
    This is the class that does the perpetual updating

        ####To build the class the following is needed:

          * **data_obj**: Object of class ATT (see above)
          * **ids_obj**: Object of class IdDict (see above)

        ####During the initializing process the following self.objects are generated:

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
          * **self.mrca_ncbi**: ncbi identifier of mrca

          * **self.tmpfi**: path to a file or folder???
          * **self.blast_subdir**: path to folder that contains the files writen during blast

          * **self.newseqs_file**: filename of files that contains the sequences from self.new_seqs_otu_id
          * **self.date**: Date of the run - may lag behind real date!
          * **self.repeat**: either 1 or 0, it is used to determine if we continue updating the tree, no new seqs found = 0
          * **self.newseqs_acc**: list of all gi_ids that were passed into remove_identical_seq(). Used to speed up adding process
          * **self.blacklist**: list of gi_id of sequences that shall not be added or need to be removed. Supplied by user.
          * **self.acc_list_mrca**: list of all gi_ids available on GenBank for a given mrca. Used to limit possible seq to add.
          * **self.seq_filter**: list of words that may occur in otu_dict.status and which shall not be used
                        in the building of FilterBlast.sp_d (that's the main function), but it is also used as assert
                        statement to make sure unwanted seqs are not added.
          * **self.unpublished**: True/False. Used to look for local unpublished seq that shall be added if True.
          * **self.path_to_local_seq:** Usually False, contains path to unpublished sequences if option is used.

        ####Following functions are called during the init-process:
        * **self.reset_markers()**: 
             adds things to self: I think they are used to make sure certain function run, if program crashed and pickle file is read in.
            * self._blasted: 0/1, if run_blast() was called, it is set to 1 for the round.
            * self._blast_read: 0/1, if read_blast_wrapper() was called, it is set to 1 for the round.
            * self._identical_removed: 0
            * self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
            * self._query_seqs_aligned: 0
            * self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
            * self._reconciled: 0
            * self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
        * **self.OToL_unmapped_tips()**: function that either removes or mappes unmapped taxa from OToL studies
    """

    def __init__(self, data_obj, ids_obj):
        # todo check input types assert()
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = self.ids.config
        self.new_seqs = {}
        self.new_seqs_otu_id = {}
        self.otu_by_gi = {}
        self._to_be_pruned = []
        self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "tmp.fasta"
        self.date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        self.repeat = 1
        self.newseqsgi = []
        self.blacklist = []
        self.gi_list_mrca = []
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,"]
        self.reset_markers()
        if self.config.blast_loc == 'local' and len(self.gi_list_mrca) == 0:
            self.gi_list_mrca = self.get_all_gi_mrca()

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

    def run_blast(self, delay = 14):  # TODO Should this be happening elsewhere?
        """generates the blast queries and saves them to xml"""
        debug("run_blast")
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            # TODO temp until I fix delete
            if otu_id in self.data.otu_dict:
                if _VERBOSE:
                    sys.stdout.write("blasting {}\n".format(otu_id))
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days)
                if time_passed > delay:  # TODO make configurable
                    equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi,
                                                                   last_blast,
                                                                   today)
                    query = seq.symbols_as_string().replace("-", "").replace("?", "")
                    if self.config.gifilename is True:
                        xml_fi = "{}/{}.xml".format(self.blast_subdir, self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label))
                    else:
                        xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
                    if _DEBUG:
                        sys.stdout.write("attempting to write {}\n".format(xml_fi))
                    if not os.path.isfile(xml_fi):
                        if _VERBOSE:
                            sys.stdout.write("blasting seq {}\n".format(taxon.label))
                        if self.config.blast_loc == 'local':
                            fi_old = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                            fi_old.write(">{}\n".format(taxon.label))
                            fi_old.write("{}\n".format(query))
                            fi_old.close()
                            blastcmd = "blastn -query " + \
                                       "{}/tmp.fas".format(self.blast_subdir) + \
                                       " -db {}nt -out ".format(self.config.blastdb) + \
                                       xml_fi + \
                                       " -outfmt 5 -num_threads {}".format(self.config.num_threads) + \
                                       " -max_target_seqs  {} -max_hsps {}".format(self.config.hitlist_size, self.config.hitlist_size) #TODO query via stdin
                            os.system(blastcmd)
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        if self.config.blast_loc == 'remote':
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
                            save_file = open(xml_fi, "w")
                            save_file.write(result_handle.read())
                            result_handle.close()
                            save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                    # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                    #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
                    else:
                        # changes date of blasted accordingly, if file is already present in the folder
                        if _DEBUG:
                            sys.stdout.write("file {} exists in current blast run. Will not blast. delete file to force\n".format(xml_fi))
                        if _DEBUG_MK == 1:
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                else:
                    if _VERBOSE:
                        sys.stdout.write("otu {} was last blasted {} days ago and is not being re-blasted."
                                          "Use run_blast(delay = 0) to force a search.\n".format(otu_id, last_blast))
        self._blasted = 1

    def get_all_gi_mrca(self):
        """get all available gi numbers from Genbank for mrca.
        The list will be used to filter out sequences from the local Blast search,
        that do not belong to ingroup."""
        # gi list limited to 100000000, for huge trees that is a problem
        debug("get_all_gi_mrca")
        Entrez.email = self.config.email
        handle = Entrez.esearch(db="nucleotide", term="txid{}[Orgn]".format(self.mrca_ncbi),
                                usehistory='n', RetMax=100000000)
        records = Entrez.read(handle)
        id_list = records['IdList']
        id_list = [int(x) for x in id_list]
        return id_list

    def read_blast(self, blast_dir=None):
        """reads in and prcesses the blast xml files"""
        debug("read blast")
        if blast_dir:
            if _VERBOSE:
                sys.stdout.write("blast dir is {}\n".format(blast_dir))
            self.blast_subdir = os.path.abspath(blast_dir)
        else:
            if _VERBOSE:
                sys.stdout.write("blast dir is {}\n".format(self.blast_subdir))
            if not os.path.exists(self.blast_subdir):
                os.mkdir(self.blast_subdir)
        if not self._blasted:
            self.run_blast()
        assert os.path.exists(self.blast_subdir)
        # because local db has no taxon info, needed to limit to group of interest
        # has been moved
        # if self.config.blast_loc == 'local' and len(self.gi_list_mrca) == 0:
        #     self.gi_list_mrca = self.get_all_gi_mrca()
            # debug("ignore mrca gi for now")
        for taxon in self.data.aln:
            # debug("add blast seq to new seqs")
            if self.config.gifilename is True:
                xml_fi = "{}/{}.xml".format(self.blast_subdir,self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label))
            else:
                xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
            if _DEBUG:
                sys.stdout.write("attempting to read {}\n".format(xml_fi))
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    if _VERBOSE:
                        sys.stdout.write(".")
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    gi_id = int(alignment.title.split('|')[1])
                                    if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                                        pass
                                    else:
                                        if gi_id not in self.data.gi_dict:  # skip ones we already have
                                            self.new_seqs[gi_id] = hsp.sbjct
                                            self.data.gi_dict[gi_id] = alignment.__dict__
                except ValueError:
                    sys.stderr.write("Problem reading {}, skipping\n".format(xml_fi))
        self.date = str(datetime.date.today())
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from GenBank after evalue filtering\n".format(len(self.new_seqs)))

        self._blast_read = 1

    # TODO this should go back in the class and should prune the tree

    def get_sp_id_of_otulabel(self, label):
        """gets the species name and the corresponding ncbi id of the otu
        """
        spn_of_label = None
        if '^ot:ottTaxonName' in self.data.otu_dict[label].keys():
            spn_of_label = self.data.otu_dict[label]['^ot:ottTaxonName']
        elif '^user:TaxonName' in self.data.otu_dict[label].keys():
            spn_of_label = self.data.otu_dict[label]['^user:TaxonName']
        if spn_of_label is not None:
            spn_of_label = str(spn_of_label).replace(" ", "_").replace("-", "_")
        if '^ncbi:taxon' in self.data.otu_dict[label].keys():
            id_of_label = self.data.otu_dict[label]['^ncbi:taxon']
        elif spn_of_label in self.ids.otu_rank.keys():
            id_of_label = int(self.ids.otu_rank[spn_of_label]["taxon id"])
        else:
            tax_name = self.ids.get_rank_info(taxon_name=spn_of_label)
            id_of_label = int(self.ids.otu_rank[tax_name]["taxon id"])
        return id_of_label

    def seq_dict_build(self, seq, label, seq_dict):  # Sequence needs to be passed in as string.
        """takes a sequence, a label (the otu_id) and a dictionary and adds the
        sequence to the dict only if it is not a subsequence of a
        sequence already in the dict.
        If the new sequence is a super sequence of one in the dict, it
        removes that sequence and replaces it
        """
        # TODO unify spp name somehow?
        id_of_label = self.get_sp_id_of_otulabel(label)
        if _DEBUG:
            sys.stdout.write("taxon is ott{}\n".format(id_of_label))
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        continue_search = False
        never_add = False

        for tax_lab in tax_list:
            existing_id = self.get_sp_id_of_otulabel(tax_lab)
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(inc_seq) >= len(new_seq):  # if seq is identical and shorter
                if inc_seq.find(new_seq) != -1:
                    # if (existing_taxa != spn_of_label and existing_taxa is not  None) or 
                    if type(existing_id) == int and existing_id != id_of_label:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, but different species name\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; subsequence, but different species"
                        seq_dict[label] = seq
                        if _DEBUG_MK == 1:
                            print(id_of_label, "and", existing_id, "subsequences, but different sp. concept")
                        continue_search = True
                        continue
                    else: # subseq of same sp.
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        if _DEBUG_MK == 1:
                            print(id_of_label, " not added, subseq of ", existing_id)
                        never_add = True
                        continue
                    return seq_dict
            else:  # if seq is longer and identical
                if new_seq.find(inc_seq) != -1:
                    if self.data.otu_dict[tax_lab].get('^physcraper:status') == "original":
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        if _DEBUG_MK == 1:
                            print(id_of_label, " and ", existing_id, "added")
                        continue_search = True
                        continue
                    elif type(existing_id) == int and existing_id != id_of_label:
                        # elif spn_of_label not in exists:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, but different species concept\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; supersequence, but different species"
                        seq_dict[label] = seq
                        if _DEBUG_MK == 1:
                            print(id_of_label, "and", existing_id, "supersequence, but different sp. concept")
                        continue_search = True
                        continue
                    else:
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added in place of {}".format(tax_lab)
                        if _DEBUG_MK == 1:
                            print(id_of_label, "added, instead of ", existing_id)
                        continue_search = True
                        continue
                    return seq_dict

        if continue_search is True or never_add is True:
            if (self.data.otu_dict[label]['^physcraper:status'].split(' ')[0] in self.seq_filter) or never_add is True:
                try:
                    del seq_dict[label]
                except:
                    debug("label was never added to seq_dict")
                try:
                    self.data.remove_taxa_aln_tre(label)
                except:
                    debug("label was never added to aln or tre")
                self.data.otu_dict[label]['^physcraper:status'] = "removed in seq dict build"  # should not be the word 'deleted', as this is used in self.seq_filter
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
        Does not test if they are identical to ones in the original alignment....
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
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        assert self.config.seq_len_perc <= 1
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        for gi, seq in self.new_seqs.items():
            if self.blacklist is not None and gi in self.blacklist:
                pass
                if _DEBUG:
                    sys.stdout.write("should not be hitting blacklist, not implemented\n")
            elif gi in self.newseqsgi:  # added to increase speed. often seq was found in another blast file
                pass
            else:
                self.newseqsgi.append(gi)
                if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                    if _DEBUG:
                        sys.stdout.write("checking gi {} for similarity\n".format(gi))
                    otu_id = self.data.add_otu(gi, self.ids)
                    self.seq_dict_build(seq, otu_id, tmp_dict)
        for tax in old_seqs:
           # try:
                del tmp_dict[tax]
           # except KeyError:
           #     pass
        self.new_seqs_otu_id = tmp_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from genbank after removing identical seq, of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def dump(self, filename=None):
        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/scrape_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi)
        # TODO... write as proper nexml?!

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
        """runs papara on the tree, the alinment and the new query sequences"""
        cwd = os.getcwd()
        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        if _VERBOSE:
            sys.stdout.write("aligning query sequences \n")
        # note: sometimes there are still sp in any of the aln/tre and I still have not found out why sometimes the label is needed
        for tax_lab in self.data.aln.taxon_namespace:
            if tax_lab not in self.data.tre.taxon_namespace:
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_lab in self.data.tre.taxon_namespace:
            if tax_lab not in self.data.aln.taxon_namespace:
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_aln in self.data.aln.taxon_namespace:
            if tax_aln not in self.data.tre.taxon_namespace:
                self.data.remove_taxa_aln_tre(tax_aln.label)
        for tax_tre in self.data.tre.taxon_namespace:
            if tax_tre not in self.data.aln.taxon_namespace:
                self.data.remove_taxa_aln_tre(tax_tre.label)
        # hack for the alien taxa thing
        self.data.reconcile()
        self.data.write_papara_files()
        os.chdir(self.workdir)  # Clean up dir moving
        try:
            assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
            subprocess.call(["papara",
                             "-t", "random_resolve.tre",
                             "-s", "aln_ott.phy",
                             "-q", self.newseqs_file,
                             "-n", papara_runname])  # FIx directory ugliness
            if _VERBOSE:
                sys.stdout.write("Papara done")
            # self.data.rewrite_files(inputfn="papara_alignment.extended")
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it installed?\n")
                sys.exit()
            # handle file not found error.
            else:
                # Something else went wrong while trying to run `wget`
                raise
        os.chdir(cwd)
        assert os.path.exists(path="{}/papara_alignment.{}".format(self.workdir, papara_runname))
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        self.data.reconcile()
        self._query_seqs_aligned = 1

    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new quesry seqs
        Just for placement, to use as starting tree."""
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

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree"""
        cwd = os.getcwd()
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        # TODO MK: work on it, first step of not using starting tree was wrong, if that is working un-comment the following stuff
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
        # run bootstrap
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "papara_alignment.extended",
                         "-p", "1", "-b", "1", "-#", "autoMRE",
                         "-n", "{}".format(self.date)])
        # make bipartition tree
        # is the -f b command
        # -z specifies file with multiple trees
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

    def remove_blacklistitem(self):
        """This removes items from aln, and tree, if the corresponding GEnbank identifer were added to the blacklist.

        Note, that seq that were not added because they were similar to the one being removed here, are lost
        (that should not be a major issue though, as in a new blast_run, they can be added.)
        """
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            if gi_id in self.blacklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, gi is part of blacklist"
        self.data.reconcile()

    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignment with the expanded ones"""

        if self.blacklist:
            self.remove_blacklistitem()
        if len(self.new_seqs) > 0:
            self.remove_identical_seqs()  # Running this twice in a reow removes all seqs before adding to alignment...
            self.data.write_files()  # should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:  # TODO rename to something more intutitive
                self.write_query_seqs()
                self.align_query_seqs()
                self.data.reconcile()
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
                if os.path.exists("{}/last_completed_update".format(self.workdir)):
                    os.rename(self.tmpfi, "{}/last_completed_update".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    # debug(filename.split("/")[-1])
                    if self.config.gifilename is not True:
                        if not os.path.exists("{}/previous_run".format(self.workdir)):
                            os.makedirs('{}/previous_run/'.format(self.workdir))
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                    else:
                        if not os.path.exists("{}/previous_run".format(self.workdir)):
                            os.makedirs('{}/previous_run/'.format(self.workdir))
                        os.rename(filename, "{}/{}".format(self.blast_subdir, filename.split("/")[-1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file),
                          "{}/previous_run/newseqs.fasta".format(self.workdir))
                try:
                    self.data.write_labelled(label='^ot:ottTaxonName')
                except:
                    self.data.write_labelled(label='user:TaxonName')
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {}  # Wipe for next run
                self.new_seqs_otu_id = {}
                self.newseqsgi = []
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
        localfiles = os.listdir(path_to_local_seq)
        for index, item in enumerate(localfiles):
            item = str(item)
            if item.startswith(".~"):
                localfiles[index] = None
        localfiles = filter(None, localfiles)
        # TODO MK: add assert to ensure that every file is a fasta file in the folder
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
                local_blast.write_blast_files(self.workdir, key, seq, db=True, fn="local_unpubl_seq")
        os.chdir(os.path.join(self.workdir, "blast"))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
        os.system(cmd1)

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
                key = species name/id
                value = dictionary:
                    key = otuID
                    value = otu_dict entry
        self.sp_seq_d: dictionary
                key = species name/id
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

    def sp_dict(self, downtorank=None):
        """Takes the information from the Physcraper otu_dict and makes a dict with species name as key and
        the corresponding seq information from aln and blast seq, it returns self.sp_d.

        This is generated to make information for the filtering class more easily available. self.sp_d sums up which
        information are available per taxonomic concept and which have not already been removed during either
        the remove_identical_seq steps or during a filtering run of an earlier cycle.
        """
        # Note: has test, test_sp_d.py, runs
        self.downtorank = downtorank
        debug("make sp_dict")
        self.sp_d = {}
        for key in self.data.otu_dict:
            if self.data.otu_dict[key]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug(self.downtorank)
                if '^ot:ottTaxonName' in self.data.otu_dict[key]:
                    spn = self.data.otu_dict[key]['^ot:ottTaxonName']
                elif '^user:TaxonName' in self.data.otu_dict[key]:
                    spn = self.data.otu_dict[key]['^user:TaxonName']
                if spn is None:
                    # debug("value is None")
                    gi_id = self.data.otu_dict[key]['^ncbi:gi']
                    # debug(gi_id)
                    # debug(type(gi_id))
                    spn = self.ids.get_rank_info(gi_id=gi_id)
                    if spn is None:
                        print("something is going wrong!Check species name")
                if self.downtorank is not None:
                    spn = str(spn).replace(" ", "_")
                    if spn not in self.ids.otu_rank.keys():
                        self.ids.get_rank_info(taxon_name=spn)
                        spn = str(spn).replace(" ", "_")
                    # if spn in self.ids.otu_rank.keys():
                    lineage2ranks = self.ids.otu_rank[spn]["rank"]
                    # else:
                    #     self.ids.get_rank_info(taxon_name=spn)
                    #     # debug(self.ids.otu_rank.keys())
                    #     lineage2ranks = self.ids.otu_rank[str(spn).replace(" ", "_")]["rank"]
                    #     # debug(lineage2ranks)s
                    ncbi = NCBITaxa()
                    for key_rank, val in lineage2ranks.iteritems():
                        if val == downtorank:
                            tax_id = key_rank
                            value_d = ncbi.get_taxid_translator([tax_id])
                            spn = value_d[int(tax_id)]
                spn = str(spn).replace(" ", "_")
                if spn in self.sp_d:
                    self.sp_d[spn].append(self.data.otu_dict[key])
                else:
                    self.sp_d[spn] = [self.data.otu_dict[key]]
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
            # loop to populate dict. key1 = sp name, key2= gb id, value = seq,
            # number of items in key2 will be filtered according to threshold and already present seq
            seq_d = {}
            tres_minimizer = 0
            for gi_id in self.sp_d[key]:
                # if statement should not be necessary as it is already filtered in the step before, I leave it in for now.
                if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    # I am using the next if to delimit which seq where already present from an earlier run,
                    # they will get a sp name (str), in order to distinguish them from newly found seq,
                    # which will have the gi (int). This differentiation is needed in the filtering blast step.
                    if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                        if '^user:TaxonName' in gi_id:
                            user_name = gi_id['^user:TaxonName']
                        elif '^ot:ottTaxonName' in gi_id:
                            user_name = gi_id['^ot:ottTaxonName']
                        for user_name_aln, seq in self.data.aln.items():
                            if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
                                otu_dict_label = self.data.otu_dict[user_name_aln.label]['^user:TaxonName']
                            elif '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
                                otu_dict_label = self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']
                            if user_name == otu_dict_label:
                                seq = seq.symbols_as_string().replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[user_name_aln.label] = seq
                        # if '^user:TaxonName' in gi_id:
                        #     # tres_minimizer += 1
                        #     user_name = gi_id['^user:TaxonName']
                        #     for user_name_aln, seq in self.data.aln.items():
                        #         if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
                        #             otu_dict_label = self.data.otu_dict[user_name_aln.label]['^user:TaxonName']
                        #             if user_name == otu_dict_label:
                        #                 seq = seq.symbols_as_string().replace("-", "")
                        #                 seq = seq.replace("?", "")
                        #                 seq_d[user_name_aln.label] = seq
                        # elif '^ot:ottTaxonName' in gi_id:
                        #     # tres_minimizer += 1
                        #     user_name = gi_id['^ot:ottTaxonName']
                        #     for user_name_aln, seq in self.data.aln.items():
                        #         if '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
                        #             otu_dict_label = self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']
                        #             if user_name == otu_dict_label:
                        #                 seq = seq.symbols_as_string().replace("-", "")
                        #                 seq = seq.replace("?", "")
                        #                 seq_d[user_name_aln.label] = seq
                    else:
                        if "^ncbi:gi" in gi_id:  # this should not be needed: all new blast seq have gi
                            # debug("gi in gi_id")
                            gi_num = int(gi_id['^ncbi:gi'])
                            if gi_num in self.new_seqs.keys():
                                seq = self.new_seqs[gi_num]
                                seq = seq.replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[gi_num] = seq
                    self.sp_seq_d[key] = seq_d
        # debug(self.sp_seq_d)
        return

    def run_local_blast(self, blast_seq, blast_db, output=None):
        """Runs  a local blast to get a measurement of differentiation between available sequences for the same taxon concept.

        The blast run will only be run if there are more sequences found than specified by the threshold value.
        When several sequences are found per taxon, blasts each seq against all other ones found for that taxon.
        The measure of differentiation will then be used to be able to select a random representative from the taxon concept,
        but allows to exclude potential mis-identifications.
        In a later step (select_seq_by_local_blast) sequences will be selected based on the blast results generated here.

        """
        # Note: has test, runs -> test_run_local_blast.py
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        out_fn = "{}_tobeblasted".format(str(blast_seq))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(blast_seq)
        os.system(cmd1)
        if output is None:
            cmd2 = "blastn -query {} -db {}_db -out output_{}.xml -outfmt 5".format(out_fn, blast_db, out_fn)
        else:
            cmd2 = "blastn -query {} -db {}_db -out {} -outfmt 5".format(out_fn, blast_db, output)
        os.system(cmd2)
        os.chdir(general_wd)

    def calculate_mean_sd(self, hsp_scores):
        """Calculates standard deviation, mean of scores which are used as a measure of sequence differentiation
        for a given taxonomic concept.

        This is being used to select a random representative of a taxonomic concept later.
        """
        # Note: has test, runs: test_calculate_mean_sd.py
        debug('calculate_mean_sd')
        total_seq = 0
        bit_sum = 0
        bit_l = []
        for gi_num in hsp_scores:
            total_seq += 1
            bit_sum += hsp_scores[gi_num]["hsp.bits"]
            bit_l.append(hsp_scores[gi_num]["hsp.bits"])
        bit_sd = float(numpy.std(bit_l))
        mean_hsp_bits = float(bit_sum / total_seq)
        mean_sd = {"mean": mean_hsp_bits, "sd": bit_sd}
        return mean_sd

    def read_local_blast(self, seq_d, fn):
        """Reads the files of the local blast run and returns sequences below a value
        (within the std of the mean scores of the hsp.bit scores at the moment).

        (this is to make sure seqs chosen are representative of the taxon)
        """
        # Note: has test, runs: test_read_local_blast.py
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        output_blast = "output_{}_tobeblasted.xml".format(fn)
        xml_file = open(output_blast)
        os.chdir(general_wd)
        blast_out = NCBIXML.parse(xml_file)
        hsp_scores = {}
        for record in blast_out:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    gi_id = alignment.title.split(" ")[1]
                    try:
                        gi_id = int(gi_id)
                    except:
                        gi_id = gi_id
                    hsp_scores[gi_id] = {"hsp.bits" : hsp.bits, "hsp.score" : hsp.score, "alignment.length" : alignment.length, "hsp.expect" : hsp.expect}
        # make values to select for blast search, calculate standard deviation,mean
        mean_sd = self.calculate_mean_sd(hsp_scores)
        # select which sequences to use
        seq_blast_score = {}
        for gi_id in hsp_scores:  # use only seq that are similar to mean plus minus sd
            if (hsp_scores[gi_id]['hsp.bits'] >= mean_sd['mean'] - mean_sd['sd']) & (hsp_scores[gi_id]['hsp.bits'] <= mean_sd['mean'] + mean_sd['sd']):
                if gi_id in seq_d:
                    seq_blast_score[gi_id] = seq_d[gi_id]
        return seq_blast_score

    def select_seq_by_local_blast(self, seq_d, fn, treshold, count):
        """Selects number of sequences from local_blast to fill up to the threshold. It returns a filtered_seq dictionary.

        It will only include species which have a blast score of mean plus/minus sd.
        Is used after read_local_blast.
        """
        # Note: has test,test_select_seq_by_local_blast.py, runs
        debug("select_seq_by_local_blast")
        seq_blast_score = self.read_local_blast(seq_d, fn)
        random_seq_ofsp = {}
        if (treshold - count) <= 0:
            debug("already too many samples of sp in aln, skip adding more.")
        elif len(seq_blast_score.keys()) == (treshold - count):
            random_seq_ofsp = seq_blast_score
        elif len(seq_blast_score.keys()) > (treshold - count):
            random_seq_ofsp = random.sample(seq_blast_score.items(), (treshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        elif len(seq_blast_score.keys()) < (treshold - count):
            random_seq_ofsp = seq_blast_score
        if len(random_seq_ofsp) > 0:
            for key, val in random_seq_ofsp.items():
                self.filtered_seq[key] = val
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
        #
        # code is untested, tested code is the code below in the comment
        seq_w_maxlen = {}
        for key, val in self.sp_seq_d[taxon_id].iteritems():
            if len(val) == len(max_len):
                seq_w_maxlen[key] = val

        if (treshold - count) <= 0:
            debug("already to many samples of sp in aln, skip adding more.")
            random_seq_ofsp = None
        elif len(seq_w_maxlen) == (treshold - count):
            random_seq_ofsp = seq_w_maxlen
        elif len(seq_w_maxlen) > (treshold - count):
            random_seq_ofsp = random.sample(seq_w_maxlen.items(), (treshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        else:
            toselect = range(len(seq_w_maxlen), (treshold - count))
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
        if random_seq_ofsp is not None:
            for key in random_seq_ofsp.keys():
                # debug(key)
                self.filtered_seq[key] = random_seq_ofsp[key]
        # print(self.filtered_seq)
        # """select new sequences by length"""
        # print("do something")
        # for key in self.sp_seq_d:
        #     count = 0
        #     if len(self.sp_seq_d[key]) > treshold:
        #         print(key)
        #         for sp_keys in self.sp_seq_d[key].keys():
        #             if isinstance(sp_keys, str) == True:
        #                 count += 1
        #         max_len = max(self.sp_seq_d[key].values())
        #         ### !!! sometimes the only seq in seq_blast_score is the
        #         original seq, then this is the one to be added, but it will be removed,
        #         # later as it is no new seq! thus no new seq for that species is added
        #         ##
        #         seq_w_maxlen = {}
        #         for k, v in self.sp_seq_d[key].iteritems():
        #             if len(v) == len(max_len):
        #                 seq_w_maxlen[k] = v
        #         if (treshold - count) <= 0:
        #             print("already to many samples of sp in aln, skip adding more.")
        #             random_seq_ofsp = None
        #         elif len(seq_w_maxlen) == (treshold - count):
        #             random_seq_ofsp = seq_w_maxlen
        #         elif len(seq_w_maxlen) > (treshold - count):
        #             random_seq_ofsp = random.sample(seq_w_maxlen.items(), (treshold - count))
        #             random_seq_ofsp = dict(random_seq_ofsp)
        #         else:
        #             toselect = range(len(seq_w_maxlen), (treshold - count))
        #             keymax = seq_w_maxlen.keys()
        #             subdict = {k:v for k, v in self.sp_seq_d[key].iteritems() if k not in keymax}
        #             2len = max(subdict.values())
        #             seq2len = {}
        #             for k, v in subdict.iteritems():
        #                 if len(v) == len(2len):
        #                     seq2len[k] = v
        #             random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
        #             random_seq_ofsp = dict(random_seq_ofsp)
        #             random_seq_ofsp.update(seq_w_maxlen)
        #         if random_seq_ofsp != None:
        #             for key in random_seq_ofsp.keys():
        #                 # print(key)
        #                 self.filtered_seq[key] = random_seq_ofsp[key]
        # # print(self.filtered_seq)

    def add_all(self, key):
        """It adds all seq to filtered_seq dict as the number of sequences present is smaller than the threshold value.

        It is only used, when sequences selection happens via blasting.

        Note: has test, test_add_all.py

        :param key: key of self.sp_d (taxon name)
        :return: self.filtered_seq
        """
        debug('add_all')
        for gi_id in self.sp_d[key]:
            if '^physcraper:status' in gi_id:
                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if gi_id['^physcraper:last_blasted'] == '1800/01/01':
                        gi_num = gi_id['^ncbi:gi']
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
        for gi_id in self.sp_d[key]:
            # this if should not be necessary
            if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                if '^user:TaxonName' in gi_id:
                    spn_name = gi_id['^user:TaxonName'].replace(" ", "_")
                elif '^ot:ottTaxonName' in gi_id:
                    spn_name = gi_id['^ot:ottTaxonName'].replace(" ", "_")
                for spn_name_aln, seq in self.data.aln.items():
                    if '^user:TaxonName' in self.data.otu_dict[spn_name_aln.label]:
                        otu_dict_name = self.data.otu_dict[spn_name_aln.label]['^user:TaxonName']
                    elif '^ot:ottTaxonName' in self.data.otu_dict[spn_name_aln.label]:
                        otu_dict_name = self.data.otu_dict[spn_name_aln.label]['^ot:ottTaxonName']
                    if spn_name == otu_dict_name:
                        nametoreturn = spn_name_aln.label
                    # for spn_name_aln in self.data.aln.keys():
                    #     if '^ot:ottTaxonName' in self.data.otu_dict[spn_name_aln.label]:
                    #         if spn_name == self.data.otu_dict[spn_name_aln.label]['^ot:ottTaxonName']:
                    #             nametoreturn = spn_name_aln.label
            # the next lines where added because the test was breaking,
            # needs thorough testing if it not breaks something else now.
            if nametoreturn is not None:
                try:
                    nametoreturn = spn_name.replace(" ", "_")
                except:
                    debug("do something?")

            else:
                break

        return nametoreturn


    def loop_for_write_blast_files(self, key):
        """This loop is needed to be able to write the local blast files for the filtering step correctly.

        Function returns a filename for the filter blast, which were generated with 'get_name_for_blastfiles'.
        """
        # Note: has test,test_loop_for_blast.py: runs
        # debug("length of sp_d key")
        # debug(len(self.sp_d[key]))
        nametoreturn = self.get_name_for_blastfiles(key)
        for gi_id in self.sp_d[key]:
            # debug("in writing file for-loop")
            # this if should not be necessary
            if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug("generate files used for blast")
                if gi_id['^physcraper:last_blasted'] != '1800/01/01':  # old seq
                    if '^user:TaxonName' in gi_id:
                        spn_name = gi_id['^user:TaxonName']
                    elif '^ot:ottTaxonName' in gi_id:
                        spn_name = gi_id['^ot:ottTaxonName']
                    for spn_name_aln, seq in self.data.aln.items():
                        if '^user:TaxonName' in self.data.otu_dict[spn_name_aln.label]:
                            otu_dict_name = self.data.otu_dict[spn_name_aln.label]['^user:TaxonName']
                        elif '^ot:ottTaxonName' in self.data.otu_dict[spn_name_aln.label]:
                            otu_dict_name = self.data.otu_dict[spn_name_aln.label]['^ot:ottTaxonName']
                        if spn_name == otu_dict_name:
                            # if selectby == "blast":
                            filename = spn_name_aln.label
                            if self.downtorank is not None:
                                filename = key
                            self.write_blast_files(nametoreturn, seq)
                    # if '^user:TaxonName' in gi_id:
                    #     spn_name = gi_id['^user:TaxonName']
                    #     for spn_name_aln, seq in self.data.aln.items():
                    #         if '^user:TaxonName' in self.data.otu_dict[spn_name_aln.label]:
                    #             if spn_name == self.data.otu_dict[spn_name_aln.label]['^user:TaxonName']:
                    #                 # if selectby == "blast":
                    #                 filename = spn_name_aln.label
                    #                 if self.downtorank is not None:
                    #                     filename = key
                    #                 self.write_blast_files(filename, seq)
                    #                 # blastfile_taxon_names[spn_name] = spn_name_aln.label
                    # elif '^ot:ottTaxonName' in gi_id:
                    #     spn_name = gi_id['^ot:ottTaxonName']
                    #     for spn_name_aln, seq in self.data.aln.items():
                    #         if '^ot:ottTaxonName' in self.data.otu_dict[spn_name_aln.label]:
                    #             if spn_name == self.data.otu_dict[spn_name_aln.label]['^ot:ottTaxonName']:
                    #                 # if selectby == "blast":  # should be obsolete now?!
                    #                 filename = spn_name_aln.label
                    #                 if self.downtorank is not None:
                    #                     filename = key
                    #                 self.write_blast_files(filename, seq)
                    #                 # blastfile_taxon_names[spn_name] = spn_name_aln.label
                else:
                    debug("make gilist as local blast database")
                    if "^ncbi:gi" in gi_id:
                        gi_num = int(gi_id['^ncbi:gi'])
                        debug(gi_num)
                        # if selectby == "blast":  # should be obsolete now?!
                        # debug("new")
                        file_present = False
                        # debug(gi_num)
                        if gi_num in self.new_seqs.keys():
                            file_present = True
                        if file_present:  # short for if file_present == True
                            if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                                # if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                                filename = gi_num
                                seq = self.sp_seq_d[key][gi_num]
                                if self.downtorank is not None:
                                    filename = key
                                    nametoreturn = key
                                self.write_blast_files(filename, seq, db=True, fn=nametoreturn)
                                # blastfile_taxon_names[gi_num] = gi_num
                    namegi = key
        if self.downtorank is not None:
            nametoreturn = key
        if nametoreturn is None:
            nametoreturn = namegi
        return nametoreturn

    def count_num_seq(self, taxon_id):
        """Counts how many sequences there are for a taxonomic concept,
        excluding sequences that have not been added during earlier cycles.

        This will be used for how_many_sp_to_keep.
        """
        # this counts the number of seq already added per taxonomic concept
        seq_present = 0
        if taxon_id in self.sp_seq_d.keys():
            for sp_keys in self.sp_seq_d[taxon_id].keys():
                if isinstance(sp_keys, str):
                    seq_present += 1
                if isinstance(sp_keys, unicode):
                    seq_present += 1
        # this determines if a taxonomic concept  is already present in the aln and how many new seq. were found
        new_taxon = True
        query_count = 0
        for item in self.sp_d[taxon_id]:
            if '^physcraper:status' in item and item['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                if item['^physcraper:last_blasted'] != '1800/01/01':
                    new_taxon = False
                if item['^physcraper:status'] == "query":
                    query_count += 1
        count_dict = {'seq_present': seq_present, 'query_count': query_count, 'new_taxon': new_taxon}
        return count_dict

    def how_many_sp_to_keep(self, treshold, selectby):
        """Uses the sp_seq_d and places the number of sequences according to threshold into the filterdseq_dict.

        This is essentially the key function of the Filter-class, it wraps up everything
        """
        debug("how_many_sp_to_keep")
        for taxon_id in self.sp_d:
            count_dict = self.count_num_seq(taxon_id)
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            if len(self.sp_d[taxon_id]) <= treshold:  # add all stuff to self.filtered_seq[gi_n]
                self.add_all(taxon_id)
            elif len(self.sp_d[taxon_id]) > treshold:  # filter number of sequences
                if taxon_id in self.sp_seq_d.keys():
                    if selectby == "length":
                        self.select_seq_by_length(self.sp_seq_d[taxon_id], treshold, seq_present)
                    elif selectby == "blast":
                        # if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
                        if 1 <= seq_present < treshold and count_dict["new_taxon"] is False and query_count != 0:
                            debug("seq_present>0")
                            if query_count + seq_present > treshold:  # species is not new in alignment, make blast with existing seq
                                taxonfn = self.loop_for_write_blast_files(taxon_id)
                                # # next loop does not seem to be used
                                # for element in self.sp_d[taxon_id]:
                                #     if '^ot:ottTaxonName' in element:
                                #         blast_seq = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                #         blast_db = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                if self.downtorank is not None:
                                    taxonfn = taxon_id

                                self.run_local_blast(taxonfn, taxonfn)
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], taxonfn, treshold, seq_present)
                            elif query_count + seq_present <= treshold:
                                self.add_all(taxon_id)
                        elif seq_present == 0 and count_dict["new_taxon"] is True and query_count >= 1:  # species is completely new in alignment
                            debug("completely new taxon to blast")
                            # species is completely new in alignment, \
                            # make blast with random species
                            for item in self.sp_d[taxon_id]:
                                if '^ncbi:gi' in item:
                                    self.data.add_otu(item['^ncbi:gi'], self.ids)
                            blast_seq = self.sp_seq_d[taxon_id].keys()[0]
                            if self.downtorank is not None:
                                str_db = taxon_id
                            else:
                                if type(blast_seq) == int:
                                    str_db = str(taxon_id)
                                else:
                                    str_db = str(blast_seq)
                            # write files for local blast first:
                            seq = self.sp_seq_d[taxon_id][blast_seq]
                            self.write_blast_files(str_db, seq)  # blast qguy
                            debug("blast db new")
                            blast_db = self.sp_seq_d[taxon_id].keys()[1:]

                            for blast_key in blast_db:
                                seq = self.sp_seq_d[taxon_id][blast_key]
                                self.write_blast_files(blast_key, seq, db=True, fn=str_db)  # local db
                            # make local blast of sequences
                            self.run_local_blast(str_db, str_db)
                            if len(self.sp_seq_d[taxon_id]) + seq_present >= treshold:
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], str_db, treshold, seq_present)
                            elif len(self.sp_seq_d[taxon_id]) + seq_present < treshold:
                                self.add_all(taxon_id)
        return

    def write_blast_files(self, file_name, seq, db=False, fn=None):
        """Writes local blast files which will be read by run_local_blast.
        """
        debug("writing files")
        if not os.path.exists("{}/blast".format(self.data.workdir)):
            os.makedirs("{}/blast/".format(self.data.workdir))
        if db:
            fnw = "{}/blast/{}_db".format(self.workdir, fn)
            fi_o = open(fnw, 'a')
        else:
            fnw = "{}/blast/{}_tobeblasted".format(self.workdir, file_name)
            fi_o = open(fnw, 'w')
        fi_o.write(">{}\n".format(file_name))
        fi_o.write("{}\n".format(seq))
        fi_o.close()

    def replace_new_seq(self):
        """Function to replace self.new_seqs and self.new_seqs_otu_id with the subset of filtered sequences.

        This is the final step in the FilterBlast class, from here it goes back to PhyScraper.

        :return: subsets of self.new_seqs and self.new_seqs_otu_id
        """
        debug("replace new seq")
        keylist = self.filtered_seq.keys()
        keylist = [x for x in keylist if type(x) == int]
        seq_not_added = self.new_seqs.keys()
        seq_not_added = [x for x in seq_not_added if type(x) == int]
        reduced_new_seqs_dic = {}
        for gi_num in seq_not_added:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'not added, there are enough seq per sp in tre'
        for gi_num in keylist:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        reduced_new_seqs_dic[key] = self.filtered_seq[gi_num]
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added, as representative of taxon'
        # might not be necessary, I was missing some gi's
        # when I was replacing the original one.
        # i leave the code here for now
        # reduced_gi_dict = {k: self.data.gi_dict[k] for k in keylist}
        # debug(reduced_gi_dict)
        # self.data.gi_dict.clear()
        # self.data.gi_dict = reduced_gi_dict # data.gi_dict seems to only have newly blasted stuff
        reduced_new_seqs = {k: self.filtered_seq[k] for k in keylist}
        with open(self.logfile, "a") as log:
            log.write("{} sequences added after filtering, of {} before filtering\n".format(len(reduced_new_seqs_dic), len(self.new_seqs_otu_id)))
        self.new_seqs = deepcopy(reduced_new_seqs)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic) # ! key is not exactly same format as before in new_seqs_otu_id
        # the dict does not seem to be used for something...
        # set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return self.new_seqs

    def write_otu_info(self, downtorank):
        """Writes a table to file with taxon names and number of representatives.
        """
        sp_d = self.sp_dict(downtorank)
        sp_info = {}
        for k in sp_d:
            sp_info[k] = len(sp_d[k])
        with open('taxon_sampling.csv', 'w') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in sp_info.items():
                writer.writerow([key, value])
        # fh = open('taxon_sampling.csv', 'w')
        # fh.write(str(sp_info) + "\n")
        # fh.close()


class Settings(object):
    """A class to store all settings for PhyScraper.
    """

    def __init__(self, seqaln, mattype, trfn, schema_trf, workdir, treshold=None,
                 selectby=None, downtorank=None, spInfoDict=None, add_local_seq=None,
                 id_to_spn_addseq_json=None, configfi=None, blacklist=None):
        """Initialize the settings."""
        self.seqaln = seqaln
        self.mattype = mattype
        self.trfn = trfn
        self.schema_trf = schema_trf
        self.workdir = workdir
        self.treshold = treshold
        self.selectby = selectby
        self.downtorank = downtorank
        self.spInfoDict = spInfoDict
        self.add_local_seq = add_local_seq
        self.id_to_spn_addseq_json = id_to_spn_addseq_json
        self.configfi = configfi
        self.blacklist = blacklist
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
