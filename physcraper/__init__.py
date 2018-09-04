#!/usr/bin/env python
"""Physcraper module"""
import sys
import re
import os
import csv
import subprocess
# import time
import datetime
# import glob
import json
# import unicodedata
import configparser
import pickle
# import inspect
import random
# import logging
import collections
from copy import deepcopy
from ete2 import NCBITaxa
# from urllib2 import URLError
import physcraper.AWSWWW as AWSWWW
import numpy
import glob
from Bio.Blast import NCBIWWW, NCBIXML
# from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import Entrez  # , SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from dendropy import Tree, \
    DnaCharacterMatrix, \
    DataSet, \
    datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life, taxomachine  # taxonomy,
from peyotl.nexson_syntax import extract_tree, \
    get_subtree_otus, \
    extract_otu_nexson, \
    PhyloSchema  # extract_tree_nexson, \
# from peyotl.api import APIWrapper
import concat  # is the local concat class
import ncbi_data_parser  # is the ncbi data parser class and associated functions

_DEBUG = 1
_DEBUG_MK = 1

_VERBOSE = 1



def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
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
debug(os.path.realpath(__file__))


class ConfigObj(object):
    """Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store."""

    def __init__(self, configfi):
        if _DEBUG:
            sys.stdout.write("Building config object\n")
        # debug(configfi)
        # debug(os.path.isfile(configfi))
        assert os.path.isfile(configfi)
        config = configparser.ConfigParser()
        config.read(configfi)
        self.e_value_thresh = config['blast']['e_value_thresh']
        assert is_number(self.e_value_thresh)
        self.hitlist_size = int(config['blast']['hitlist_size'])
        self.seq_len_perc = float(config['physcraper']['seq_len_perc'])
        assert 0 < self.seq_len_perc < 1
        self.get_ncbi_taxonomy = config['taxonomy']['get_ncbi_taxonomy']
        assert os.path.isfile(self.get_ncbi_taxonomy)
        self.ncbi_dmp = config['taxonomy']['ncbi_dmp']
        # gi to taxid (according to GenBank it's not updated since 2016, even though the files seems to be newwe)
        if not os.path.isfile(self.ncbi_dmp):
            os.system("rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz {}.gz".format(self.ncbi_dmp))
            os.system("gunzip taxonomy/gi_taxid_nucl.dmp.gz")
            self.ncbi_dmp = "taxonomy/gi_taxid_nucl.dmp.gz"
        # # acc to taxid
        # self.acc2taxid = config['taxonomy']['acc2taxid']
        # if not os.path.isfile(self.acc2taxid):
        #     os.system("rsync -av ftp.ncbi.nlm.nih.gov::pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz")  # TODO: command not the right one, downloaded by hand for now
        #     os.system("gunzip taxonomy/nucl_gb.accession2taxid.gz")
        #     self.acc2taxid = "taxonomy/nucl_gb.accession2taxid"
        # # rank of lineages
        # self.rankedlineages = config['taxonomy']['rankedlineages']
        # if not os.path.isfile(self.rankedlineages):
        #     os.system("rsync -av ftp.ncbi.nlm.nih.gov::pub/taxonomy/new_taxdump/new_taxdump.tar.gz") # TODO: command not the right one, downloaded by hand for now
        #     os.system("gunzip taxonomy/rankedlineage.dmp.gz new_taxdump.tar/rankedlineage.dmp ")
        #     self.rankedlineages = "taxonomy/rankedlineages.dmp"
        self.phylesystem_loc = config['phylesystem']['location']
        assert (self.phylesystem_loc in ['local', 'api'])
        self.ott_ncbi = config['taxonomy']['ott_ncbi']
        assert os.path.isfile(self.ott_ncbi)
        self.id_pickle = os.path.abspath(config['taxonomy']['id_pickle'])  # TODO what is theis doing again?
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

        if _DEBUG:
            sys.stdout.write("{}\n".format(self.email))
            if self.blast_loc == 'remote':
                sys.stdout.write("url base = {}\n".format(self.url_base))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == 'local':
                sys.stdout.write("local blast db {}\n".format(self.blastdb))


# ATT is a dumb acronym for Alignment Tree Taxa object
def get_dataset_from_treebase(study_id,
                              phylesystem_loc='api'):
    try:
        nexson = get_nexson(study_id, phylesystem_loc)
    except:
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
                                  phylesystem_loc='api'):
    """gathers together tree, alignment, and study info - forces names to otu_ids.
    Outputs AlignTreeTax object.
    an alignemnt, a
    Input can be either a study ID and tree ID from OpenTree
    Alignemnt need to be a Dendropy DNA character matrix!"""
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
            sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignement. "
                             "This may indicate a mismatch between tree and alignement\n".format(tax.label))
    # need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    workdir = os.path.abspath(workdir)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)
    # newick should be bare, but alignement should be DNACharacterMatrix


def convert(data):
    """convert json 2.7 as string problem"""
    if isinstance(data, basestring):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert, data))
    else:
        return data


def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            treefile,
                            otu_json,
                            schema_trf,
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem.
    If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.
    otu_json should encode the taxon names for each tip"""
    # Note: has test -> owndata.py

    # replace ? in seqaln with - : papara handles them as different characters
    with open(seqaln, 'r') as file:
        filedata = file.read()
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
    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        ott_ids = [otu_dict[otu].get(u'^ot:ottId', ) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)


def standardize_label(item):
    """try to make names unicode
    """
    # Note: has test, runs -> test_edit_dict_key.py
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
    except:
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
    reads input file into the var spInfo, tranaltes using an IdDict object
    using web to call Open tree, then ncbi if not found"""
    sys.stdout.write("Set up OtuJsonDict \n")
    sp_info_dict = {}
    with open(id_to_spn, mode='r') as infile:
        for lin in infile:
            tipname, species = lin.strip().split(',')
            # debug(tipname)
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
                # debug(name2taxid)
                if len(name2taxid.items()) >= 1:
                    if _DEBUG:
                        sys.stdout.write("found taxon {} in ncbi".format(spn))
                    ncbiid = name2taxid.items()[0][1][0]
                    ottid = id_dict.ncbi_to_ott[ncbiid]
                    ottname = id_dict.ott_to_name[ottid]
                else:
                    sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                     "Proceeding without taxon info\n".format(spn))
                    ottid, ottname, ncbiid = None, None, None
            sp_info_dict[otu_id] = {'^ncbi:taxon': ncbiid, '^ot:ottTaxonName': ottname, '^ot:ottId': ottid,
                                    '^ot:originalLabel': tipname, '^user:TaxonName': species,
                                    '^physcraper:status': 'original', '^physcraper:last_blasted': "1900/01/01"}
    return sp_info_dict


class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
    Hypothetically, all teh keys in the  otu_dict shopuld be cealn
    """

    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None, taxon_namespace=None):
        # TODO add assertions that inputs are correct type!!!
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
        self.tre = Tree.get(data=newick,
                            schema="newick",
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        assert isinstance(otu_dict, dict)
        self.otu_dict = otu_dict
        self.ps_otu = 1  # iterator for new otu IDs
        self._reconcile_names()
        self.workdir = os.path.abspath(workdir)  # TODO - is this where the workdir should live?
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca)
        self.ott_mrca = ingroup_mrca  # TODO: we only use .ott_mrca to infer mrca_ncbi. Why not using the ncbi one directly?
        self.orig_seqlen = []  # FIXME
        self.gi_dict = {}  # has all info about new blast seq  TODO: this should be part of physcraper class, as it has all blast information. Blast is not part of this class.
        self.orig_aln = alignment
        self.orig_newick = newick
        self._reconciled = False
        self.local_otu_json = None

    def _reconcile_names(self):
        """This checks that the tree "original labels" from phylsystem
        align with those found in the alignment. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!"""
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        aln_tax = set([tax for tax in self.aln.taxon_namespace])
        prune = treed_taxa ^ aln_tax
        missing = [i.label for i in prune]
        if missing:
            errmf = 'NAME RECONCILIATION Some of the taxa in the tree are not in the alignment or vice versa and will be pruned. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            sys.stderr.write(errm)
        self.aln.remove_sequences(prune)
        self.tre.prune_taxa(prune)
        for tax in prune:
            if tax.label in self.otu_dict:
                self.otu_dict[tax.label]['^physcraper:status'] = "deleted in name reconciliation"
            else:
                sys.stderr.write("lost taxon {} in name reconcilliation".format(tax.label))
            self.aln.taxon_namespace.remove_taxon(tax)
        assert (self.aln.taxon_namespace == self.tre.taxon_namespace)
        # reverse_otu_dict = {}  # seems not to be used
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
                    if self.otu_dict[otu].get('^ot:originalLabel') == tax.label or self.otu_dict[otu].get('^ot:originalLabel') == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))
                # assert tax.label in self.otu_dict
    # TODO - make sure all taxon labels are unique OTU ids.

    def prune_short(self, min_seqlen=0):
        """Sometimes in the de-concatenating of the original alignment
        taxa with no sequence are generated.
        This gets rid of those from both the tre and the alignement. MUTATOR
        """
        # has test, runs: test_prune_short.py
        debug("prune short")
        prune = []
        for tax, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) <= min_seqlen:
                prune.append(tax)
        if prune:
            # debug(prune)
            self.aln.remove_sequences(prune)
            self.tre.prune_taxa(prune)
            self.tre.prune_taxa_with_labels(prune)  # sometimes it does not delete it with the statement before. Tried to figure out why, have no clue yet.
            # self.aln.taxon_namespace.remove_taxon_label(tax)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short step due to sequence shorter than {}\n".format(min_seqlen))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in prune short"
            self.aln.taxon_namespace.remove_taxon_label(tax.label)  # raises error if not found, instead of remove_taxon
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        self.reconcile()

    def reconcile(self, seq_len_perc=0.75):
        """all missing data seqs are sneaking in, but from where?!"""
        # not only used in the beginning...is used to remove sequences that are shorter than 75%
        # assert self.aln.taxon_namespace == self.tre.taxon_namespace
        # debug("reconcile")
        prune = []
        # debug(self.orig_seqlen)
        avg_seqlen = sum(self.orig_seqlen) / len(self.orig_seqlen)
        # debug(avg_seqlen)
        seq_len_cutoff = avg_seqlen * seq_len_perc
        # debug(seq_len_cutoff)
        for tax, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(tax)
        if prune:
            # debug(prune)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and aln in reconcilation step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict.get(tax.label).get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            # debug(tax)
            # debug(tax.label)
            # debug(self.otu_dict[tax.label])
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in reconcile"
            # TODO: line above is unnecessary? get's overwritten in next line by remove_taxa_aln_tre
            self.remove_taxa_aln_tre(tax.label)
        aln_ids = set()
        for tax in self.aln:
            aln_ids.add(tax.label)

        # debug(len(self.otu_dict.keys()))
        debug(len(aln_ids))
        debug([item for item in self.otu_dict.keys() if item not in aln_ids])
        debug([item for item in aln_ids if item not in self.otu_dict.keys()])

        assert aln_ids.issubset(self.otu_dict.keys())
        treed_taxa = set()
        orphaned_leafs = set()
        # assert self.aln.taxon_namespace == self.tre.taxon_namespace
        # here leaf_nodes have taxa that were dropped before. Why do we have this anyways?
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
            if leaf.taxon.label not in aln_ids:
                self.otu_dict[leaf.taxon.label]['^physcraper:status'] = "deleted due to presence in tree but not aln ?!"
                orphaned_leafs.add(leaf)
                # TODO figure out why sometimes one of them works and not the other and vice versa
                self.tre.prune_taxa([leaf])
                # self.tre.prune_taxa_with_labels([leaf.taxon.label])
                self.tre.prune_taxa_with_labels([leaf.taxon.label])
                self.tre.prune_taxa_with_labels([leaf])
                treed_taxa.remove(leaf.taxon.label)
                # debug(self.otu_taxonlabel_problem.keys())
            # else:
            #     treed_taxa.add(leaf.taxon.label)
        # debug('treed_taxa')
        # debug(treed_taxa)
        # debug('aln_ids')
        # debug(aln_ids)
        # debug([item for item in treed_taxa if item not in aln_ids])
        # #debug([item for item in aln_ids if item not in treed_taxa])
        # debug(self.tre.taxon_namespace) # otu is gone from namespace, but in treed
        # debug(self.tre.as_string(schema="newick")) # otu is in tre, thus  not removed in remove_taxa_aln_tre
        assert treed_taxa.issubset(aln_ids)
        # for key in  self.otu_dict.keys():
        #      if key not in aln_ids:
        #           sys.stderr.write("{} was in otu dict but not alignment. it should be in new seqs...\n".format(key)
        self.trim()
        self._reconciled = 1

    def trim(self, taxon_missingness=0.75):
        """cuts off ends of alignment, maintaining similar to original seq len
        Important bc other while whole chromosomes get dragged in!"""
        debug('in trim')
        seqlen = len(self.aln[0])
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
        for i in range(seqlen, 0, -1):  # seqlen-1 cuts off last character of aln, I changed it.
            counts = {'?': 0, '-': 0}
            for tax in self.aln:
                call = self.aln[tax][i - 1].label  # changing seqlen-1 to seqlen requires that we have here i-1
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                stop = i
                break
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        aln_ids = set()
        for tax in self.aln:
            aln_ids.add(tax.label)
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
            sys.stdout.write("trimmed alignment ends to < {} missing taxa, start {}, stop {}\n".format(taxon_missingness, start, stop))
        return


    def add_otu(self, gi_id, ids_obj):
        """generates an otu_id for new sequences and adds them into the otu_dict.
        Needs to be passed an IdDict to do the mapping"""
     
        # debug("add_otu function")
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        # debug(gi_id)
        ncbi_id = None
        ott = None
        if type(gi_id) == int:
            # debug("gi_id is int")
            if gi_id in self.gi_dict.keys() and 'staxids' in self.gi_dict[gi_id].keys():
                spn = self.gi_dict[gi_id]['sscinames']
                ncbi_id = self.gi_dict[gi_id]['staxids']
                # debug(ncbi_id)
            else:
                spn = ids_obj.find_name(gi=gi_id)
                if spn is None:
                    sys.stderr.write("no species name returned for {}".format(gi_id))
                ncbi_id = ids_obj.map_gi_ncbi(gi_id)

        elif gi_id[:6] == "unpubl":
            # debug(self.gi_dict.keys())
            # debug(self.gi_dict[gi_id].keys())
            spn = self.gi_dict[gi_id]['^ot:ottTaxonName']
            ncbi_id = self.gi_dict[gi_id]['^ncbi:taxon']
            ott = self.gi_dict[gi_id]['^ot:ottId']
        else:
            sys.stderr.write("Something is wrong, I cannot add a new seq which has no gi id or is not unpublished.")


        if ncbi_id is None:
            print("ncbi is none")
            ncbi_id = ids_obj.ncbi_parser.get_id_from_name(spn)
            if type(gi_id) == int:
                print("add id to self")
                self.gi_ncbi_dict[gi_id] = tax_id
            self.ncbiid_to_spn[ncbi_id] = spn
            # deprecated for ncbi_parser
            # ncbi_id = ids_obj.otu_rank[spn]["taxon id"]
            # debug(ncbi_id)
        if ncbi_id in ids_obj.ncbi_to_ott.keys():
            # ncbi_id = int(ids_obj.map_gi_ncbi(gi_id))
            # try:
            ott = int(ids_obj.ncbi_to_ott[ncbi_id])
            # except:
        if ott is None:
            ott = "OTT_{}".format(self.ps_otu)
            self.ps_otu += 1
            # spn = str(ids_obj.ott_to_name[ott]).replace(" ", "_")  # seems to be unused

        if otu_id in self.otu_dict.keys():
            ott_name = ids_obj.ott_to_name.get(ott)
        else:
            ott_name = spn
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
        self.otu_dict[otu_id]['^ot:ottId'] = ott
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = ott_name
        # last_blasted date infos: 1800 = never blasted; 1900 = blasted 1x, not added; this century = blasted and added
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"

        if type(gi_id) != int:
            print(gi_id)
            # key = "unpubl_{}".format(gi_id)
            self.otu_dict[otu_id]['^user:TaxonName'] = self.gi_dict[gi_id]['localID']
            self.otu_dict[otu_id]['^physcraper:status'] = "local seq"

        if _DEBUG >= 2:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi_id, otu_id))
        # debug(otu_id)
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

    def write_labelled(self, label, treepath=None, alnpath=None, norepeats=True, gi_id=False):
        """output tree and alignement with human readable labels
        Jumps through abunch of hoops to make labels unique.
        NOT MEMORY EFFICIENT AT ALL"""
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
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"), self.otu_dict[taxon.label].get('^ot:ottTaxonName', "unk"))
            new_label = str(new_label).replace(' ', '_')
            if gi_id:
                gi_id = self.otu_dict[taxon.label].get('^ncbi:gi')
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
        """Writes out OTU dict as json"""
        assert schema in ['table', 'json']
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln, tre and otu_dict,
        takes a single taxon_label as input.
        """
        # note: has test, test_remove_taxa_aln_tre.py, runs, passes
        # debug('remove_taxa_aln_tre')
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
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, but it wasn't in the alignment..."

    def dump(self, filename=None):
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
    tree. The blast search later is limited to descendents of this
    mrca according to the ncbi taxonomy"""
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
        sys.stderr.write('No sampled taxa were found in the current sysnthetic tree. '
                         'Please find and input and approppriate OTT id as ingroup mrca in generate_ATT_from_files')
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
                         'approppriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit()
    return tax_id


def get_ott_ids_from_otu_dict(otu_dict):  # TODO put into data obj?
    """Get the ott ids from an otu dict object"""
    # never used
    ott_ids = []
    for otu in otu_dict:
        try:
            ott_ids.append(otu['^ot:ottId'])
        except KeyError:
            pass


#####################################

class IdDicts(object):
    """Wraps up the annoying conversions
    """
    # TODO - could - should be shared acrosss runs?! .... nooo.
    def __init__(self, config_obj, workdir):
        """Generates a series of name disambiguation dicts"""
        self.workdir = workdir
        self.config = config_obj
        assert self.config.email
        self.ott_to_ncbi = {}  # currently only used to find mcra ncbi id from mrca ott id
        self.ncbi_to_ott = {}  # used to get ott id for new Genbank query taxa
        self.ott_to_name = {}
        # self.otu_rank = {}  # deprecated for ncbi_parser
        self.gi_ncbi_dict = {}  # file id_map is not existing, is only filled by get_rank_info/ncbi_parser. is a smaller version of self.otu_rank.
        self.spn_to_ncbiid = {}  # spn to ncbiid, TODO: change as they are being fed by the ncbi_data_parser infos only
        self.ncbiid_to_spn = {}
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
        if os.path.isfile("{}/id_map.txt".format(workdir)):  # todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]
        # ncbi parser contains information about spn, tax_id, and ranks
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                                   nodes_file=self.config.ncbi_parser_nodes_fn)

    def get_ncbiid_from_tax_name(self, tax_name):
        """Get the ncbi id from the species name using ncbi web query.
        """
        if tax_name in self.spn_to_ncbiid:
            tax_id = self.spn_to_ncbiid[tax_name]
        else:
            try:
                # debug("try2")
                tries = 10
                for i in range(tries):
                    try:
                        tax_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
                        tax_id = int(tax_id)
                    except:
                        # debug("except esearch/read")
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            # debug("im going to raise")
                            raise
                    break
                # debug(tax_id)
                # debug(type(tax_id))
            except:
                # debug("except")
                ncbi = NCBITaxa()
                tax_info = ncbi.get_name_translator([tax_name])
                # debug(tax_info)
                if tax_info == {}:
                    debug("Taxon name does not match any name in ncbi. Check that name is written correctly!")
                tax_id = int(tax_info.items()[0][1][0])
        assert type(tax_id) is int
        self.spn_to_ncbiid[tax_name] = tax_id
        return tax_id

    def find_name(self, sp_dict=None, gi=None):
        """ Find the name in the sp_dict or of a gi. If not already known if will ask ncbi using the gi number.
        """
        # debug("find_name")
        spn = None
        if sp_dict:

            if '^user:TaxonName' in sp_dict:
                spn = sp_dict['^user:TaxonName']
            elif '^ot:ottTaxonName' in sp_dict:
                spn = sp_dict['^ot:ottTaxonName']
        if spn is None:
            if gi:
                gi_id = gi
            elif '^ncbi:gi' in sp_dict:
                # 
                gi_id = sp_dict['^ncbi:gi']
            else:
                sys.stderr.write("There is no name supplied and no gi available. This should not happen! Check name!")

            if gi_id in self.gi_ncbi_dict:
                ncbi_taxonid = self.gi_ncbi_dict[gi]
                
                # # deprecated for ncbi_parser
                # for key, value in self.spn_to_ncbiid.values():
                #     if value == ncbi_taxonid:
                #         spn = key
                if ncbi_taxonid in self.ncbiid_to_spn.keys:
                    spn = self.ncbiid_to_spn[ncbi_taxonid]
                else:
                    spn = self.ncbi_parser.get_name_from_id(ncbi_taxonid)
                    self.ncbiid_to_spn[ncbi_taxonid] = spn
            else:
                # debug(gi_id)
                # debug(type(gi_id))
                tries = 10
                Entrez.email = self.config.email
                for i in range(tries):
                    try:
                        # debug("find name efetch")
                        handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
                    except:
                        # debug("except efetch")
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            # debug("im going to raise")
                            raise
                    break
                read_handle = Entrez.read(handle)
                handle.close()
                spn = get_ncbi_tax_name(read_handle)
                ncbi_taxonid = get_ncbi_tax_id(read_handle)
                if sp_dict:
                    sp_dict['^ot:ottTaxonName'] = spn

                    sp_dict['^ncbi:taxon'] = ncbi_taxonid
        assert spn is not None
        spn = spn.replace(" ", "_")
        return spn

    def map_gi_ncbi(self, gi_id):
        """get the ncbi taxon id's for a gi input.
        """
        # debug("map_gi_ncbi")
        if _DEBUG == 2:
            sys.stderr.write("mapping gi {}\n".format(gi_id))
        if gi_id in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi_id])
        else:
            # debug(gi)
            # tax_name = self.get_rank_info(gi_id=gi)
            tax_name = self.find_name(gi=gi_id)
            # tax_id = self.get_ncbiid_from_tax_name(tax_name)  # uses internet query, next line, dmp file
            tax_id = self.ncbi_parser.get_id_from_name(tax_name)
            self.ncbiid_to_spn[tax_id] = tax_name
            # if type(gi_id) == int:
            #     print("map_gi: add id to self")
            #     self.gi_ncbi_dict[gi_id] = tax_id
            # debug(tax_name)
            # tax_id = self.otu_rank[tax_name]["taxon id"]
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
    """This is the class that does the perpetual updating"""

    def __init__(self, data_obj, ids_obj):
        # todo check input types assert()
        debug("start base class init")
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = self.ids.config  # this is already part of .ids, or not? Information are doubled.
        self.new_seqs = {}  # all new seq after read_blast
        self.new_seqs_otu_id = {}  # only new seq which passed remove_identical
        self.otu_by_gi = {}  # TODO: What was this intended for?
        self._to_be_pruned = []  # TODO: where do we use it?
        self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "tmp.fasta"
        self.date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        self.repeat = 1
        self.newseqsgi = []  # all ever added gi during any PhyScraper run,
        self.blacklist = []  # remove sequences by default
        self.gi_list_mrca = []
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,", "local"]
        self.reset_markers()
        if self.config.blast_loc == 'local' and len(self.gi_list_mrca) == 0:
            self.gi_list_mrca = self.get_all_gi_mrca()
            # debug(self.gi_list_mrca)
        self.unpublished = False
        self.path_to_local_seq = False

        # self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn, nodes_file=self.config.ncbi_parser_nodes_fn)
        # self.query_dict = {}  # for local blast txt files, equivalent to gi_dict.

    # TODO is this the right place for this?
    def reset_markers(self):
        self._blasted = 0
        self._blast_read = 0
        self._identical_removed = 0
        self._query_seqs_written = 0
        self._query_seqs_aligned = 0  # TODO: Where do we use it?
        self._query_seqs_placed = 0
        self._reconciled = 0  # TODO: Where do we use it?
        self._full_tree_est = 0

    def run_blast(self, delay=14):  # TODO Should this be happening elsewhere?
        """generates the blast queries and saves them to xml"""
        debug("run_blast")
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            # print(taxon, seq)
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
                    fi_old = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                    fi_old.write(">{}\n".format(taxon.label))
                    fi_old.write("{}\n".format(query))
                    fi_old.close()
                    blast_db = "local_unpubl_seq_db"
                    output = "tst_fn"
                    blastcmd = "blastn -query {}/tmp.fas -db {} -out output_{}.xml -outfmt 5".format(self.blast_subdir, blast_db, output)
                    os.system(blastcmd)
                    # self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                else:
                    if time_passed > delay:
                        equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi,
                                                                       last_blast,
                                                                       today)
                        if self.config.blast_loc == 'local':
                            file_ending = "txt"
                        else:
                            file_ending = "xml"
                        if self.config.gifilename is True:
                            fn = self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label)
                            fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                        else:
                            fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                        if _DEBUG:
                            sys.stdout.write("attempting to write {}\n".format(fn_path))
                        if not os.path.isfile(fn_path):
                            if _VERBOSE:
                                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                            if self.config.blast_loc == 'local':
                                cwd = os.getcwd()
                                os.chdir(self.config.blastdb)
                                fi_old = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                                fi_old.write(">{}\n".format(taxon.label))
                                fi_old.write("{}\n".format(query))
                                fi_old.close()
                                # this formats allows to get the taxonomic information at the same time
                                outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq stitle'"
                                # outfmt = " -outfmt 5"  # format for xml file type
                                blastcmd = "blastn -query " + \
                                           "{}/tmp.fas".format(self.blast_subdir) + \
                                           " -db {}nt -out ".format(self.config.blastdb) + \
                                           fn_path + \
                                           " {} -num_threads {}".format(outfmt, self.config.num_threads) + \
                                           " -max_target_seqs {} -max_hsps {}".format(self.config.hitlist_size, self.config.hitlist_size)  # TODO query via stdin
                                # debug(blastcmd)
                                os.system(blastcmd)
                                os.chdir(cwd)
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
                                    # debug(result_handle.read())
                                save_file = open(fn_path, "w")
                                save_file.write(result_handle.read())
                                result_handle.close()
                                save_file.close()
                                self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                        #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
                        else:
                            # changes date of blasted accordingly, if file is already present in the folder
                            if _DEBUG:
                                sys.stdout.write(
                                    "file {} exists in current blast run. Will not blast, delete file to force\n".format(fn_path))
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
        """reads in and processes the blast xml files"""
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
            # xml_fi = "{}/blast/{}.xml".format(self.workdir, "output_tst_fn")
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
                                fake_gi = "unpubl_{}".format(gi_id)
                                # self.make_otu_dict_entry_unpubl(fake_gi)
                                self.new_seqs[fake_gi] = hsp.sbjct
                                # debug(gi_id)
                                # debug(self.data.local_otu_json['otu{}'.format(gi_id)])
                                self.data.gi_dict[fake_gi] = {'accession': "000000{}".format(gi_counter), 'title': "unpublished", 'localID': gi_id}
                                self.data.gi_dict[fake_gi].update(self.data.local_otu_json['otu{}'.format(gi_id)])
                                gi_counter += 1
                                # self.data.gi_dict[fake_gi] = alignment.__dict__
        else:
            if not self._blasted:
                self.run_blast()
            assert os.path.exists(self.blast_subdir)
            for taxon in self.data.aln:
                # debug(taxon)
                # debug("add blast seq to new seqs")
                if self.config.blast_loc == 'local':
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                if self.config.gifilename is True:
                    fn = self.data.otu_dict[taxon.label].get('^ncbi:gi', taxon.label)
                    fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                else:
                    fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                if _DEBUG:
                    sys.stdout.write("attempting to read {}\n".format(fn_path))
                if os.path.isfile(fn_path):
                    if self.config.blast_loc == 'local':  # new method to read in txt format
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
                                # debug(type(gi_id))
                                debug(gi_id)
                                if len(self.gi_list_mrca) >= 1 and (gi_id not in self.gi_list_mrca):
                                    # debug("pass")
                                    pass
                                else:
                                    # debug("try to add to new seqs")
                                    if gi_id not in self.data.gi_dict:  # skip ones we already have            
                                        # debug("added")
                                        self.new_seqs[gi_id] = query_dict[key]['sseq']
                                        self.data.gi_dict[gi_id] = query_dict[key]
                    else:
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
        self.date = str(datetime.date.today())
        debug("len new seqs dict after evalue filter")
        debug(len(self.new_seqs))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from GenBank after evalue filtering\n".format(len(self.new_seqs)))

        self._blast_read = 1

    # TODO this should go back in the class and should prune the tree

    def get_sp_id_of_otulabel(self, label):
        """gets the species name and the corresponding ncbi id of the otu
        """
        # debug("get_spn_id_of_otulabel")
        # debug(label)
        # debug(self.data.otu_dict[label].keys())
        spn_of_label = self.ids.find_name(sp_dict=self.data.otu_dict[label])
        # if '^ot:ottTaxonName' in self.data.otu_dict[label].keys():
        #     spn_of_label = self.data.otu_dict[label]['^ot:ottTaxonName']
        # elif '^user:TaxonName' in self.data.otu_dict[label].keys():
        #     spn_of_label = self.data.otu_dict[label]['^user:TaxonName']
        # else:
        #     spn_of_label = None
        if spn_of_label is not None:
            # spn_of_label = str(spn_of_label).replace(" ", "_").replace("-", "_")
            spn_of_label = str(spn_of_label).replace(" ", "_")

        else:
            debug("Problem, no sp name found!")
        debug(spn_of_label)
        
        if spn_of_label in self.ids.spn_to_ncbiid:
            id_of_label = self.ids.spn_to_ncbiid[spn_of_label]
        else:
            id_of_label = self.ids.ncbi_parser.get_id_from_name(spn_of_label)
            self.ids.spn_to_ncbiid[spn_of_label] = id_of_label

        # following code should not be necessary anymore since intro of ncbi_parser 
        # if spn_of_label in self.ids.otu_rank.keys():
        #     # debug("i already have the rank info")
        #     id_of_label = int(self.ids.otu_rank[spn_of_label]["taxon id"])

        # elif '^ncbi:taxon' in self.data.otu_dict[label].keys():
        #     # debug(self.data.otu_dict[label])
        #     # print("no rank info yet")
        #     if self.data.otu_dict[label]['^ncbi:taxon'] is not None:
        #         # debug("taxon name is not none")
        #         # debug(self.data.otu_dict[label])
        #         id_of_label = int(self.data.otu_dict[label]['^ncbi:taxon'])
        #     else:
        #         # print("i have nothing")
        #         # tax_name = self.ids.get_rank_info(taxon_name=spn_of_label)
        #         id_of_label = self.ids.ncbi_parser.get_id_from_name(spn_of_label)
        #         # id_of_label = int(self.ids.otu_rank[tax_name]["taxon id"])
        #     # debug(id_of_label)
        #     # debug(type(id_of_label))
        # else:
        #     # debug(spn_of_label)
        #     debug("i have nothing to use and im in else")
        #     # tax_name = self.ids.get_rank_info(taxon_name=spn_of_label)
        #     id_of_label = self.ids.ncbi_parser.get_id_from_name(spn_of_label)
        #     # id_of_label = int(self.ids.otu_rank[tax_name]["taxon id"])
        return id_of_label

    def seq_dict_build(self, seq, label, seq_dict):  # Sequence needs to be passed in as string.
        """takes a sequence, a label (the otu_id) and a dictionary and adds the
        sequence to the dict only if it is not a subsequence of a
        sequence already in the dict.
        If the new sequence is a super sequence of one in the dict, it
        removes that sequence and replaces it
        """
        # TODO unify spp name somehow?
        # debug("seq_dict_build")

        # debug(label)
        id_of_label = self.get_sp_id_of_otulabel(label)
        # debug("id of label is{}".format(id_of_label))
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
            # print(len(inc_seq), len(new_seq))
            # debug(sum(self.data.orig_seqlen) / len(self.data.orig_seqlen))
            # print(id_of_label, existing_id)
            # if id_of_label == "irrelevant_sequence":
            #     print("do not add to aln!")
            #     debug(some)
            if len(new_seq) >= sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * 2.5:
                debug("seq not added because it is to long...")
            elif len(inc_seq) >= len(new_seq):  # if seq is identical and shorter
                if inc_seq.find(new_seq) != -1:
                    # if (existing_taxa != spn_of_label and existing_taxa is not  None) or 
                    # print(type(existing_id))
                    if type(existing_id) == int and existing_id != id_of_label:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, but different species name\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; subsequences, but different species"
                        seq_dict[label] = seq
                        debug("{} and {} are subsequences, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:  # subseq of same sp.
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        debug("{} not added, subseq of".format(existing_id))

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
                        debug("{} and  {} added".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    elif type(existing_id) == int and existing_id != id_of_label:
                        # elif spn_of_label not in exists:
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, but different species concept\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; supersequence, but different species"
                        seq_dict[label] = seq
                        debug("{} and  {} supersequence, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added in place of {}".format(tax_lab)
                        debug("{} added, instead of  {}".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    return seq_dict

        if continue_search is True or never_add is True:
            if (self.data.otu_dict[label]['^physcraper:status'].split(' ')[0] in self.seq_filter) or never_add is True:
                if label in seq_dict.keys():
                    del seq_dict[label]
                else:
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
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        # Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        assert self.config.seq_len_perc <= 1
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        # debug(self.config.seq_len_perc)
        # debug("seq cutoff")
        # debug(seq_len_cutoff)
        for gi, seq in self.new_seqs.items():
            debug(gi)
            if self.blacklist is not None and gi in self.blacklist:
                debug("gi in blacklist, not added")
                pass
            elif gi in self.newseqsgi:  # added to increase speed. often seq was found in another blast file
                debug("passed, was already added")
                pass
            else:
                debug("add to aln if not similar exists")
                if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                    
                    if type(gi) == int or gi.isdigit():
                        debug("gi is digit")
                        if type(gi) != int:
                            sys.stdout.write("WARNING: gi {} is no integer. Will convert value to int\n".format(gi))
                            debug("WARNING: gi {} is no integer. Will convert value to int\n".format(gi))
                            gi = int(gi)
                    # if self.config.blast_loc == 'local':
                    #     localblast = True
                    # else:
                    #     localblast = False
                    self.newseqsgi.append(gi)
                    otu_id = self.data.add_otu(gi, self.ids)
                    # debug(otu_id)
                    # debug("go to seq_dict_build")
                    self.seq_dict_build(seq, otu_id, tmp_dict)
            # else:
            #     debug("gi was already compared")
        for tax in old_seqs:
            try:
                del tmp_dict[tax]
            except:
                pass
        self.new_seqs_otu_id = tmp_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        debug("len new seqs dict after remove identical")
        debug(len(self.new_seqs_otu_id))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def find_otudict_gi(self):
        ncbigi_list =[]
        for key, val in self.data.otu_dict.items():
            if '^ncbi:gi' in val:
                gi_otu_dict = val["^ncbi:gi"]
                ncbigi_list.append(gi_otu_dict)
        return ncbigi_list

    def dump(self, filename=None):
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
        """runs papara on the tree, the alinment and the new query sequences"""
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
        # I still have not found out why sometimes the label is needed
        self.remove_alien_aln_tre()
        # hack for the alien taxa thing
        debug("added another reconcile step, maybe that helps with the S.glaber/alien taxa problem")
        self.data.reconcile()
        self.data.write_papara_files()
        os.chdir(self.workdir)  # Clean up dir moving
        try:
            debug("I call papara")
            assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
            # debug(self.newseqs_file)
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
        # debug(self.data.aln.taxon_namespace)
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        self.data.reconcile()
        self._query_seqs_aligned = 1

    def remove_alien_aln_tre(self):
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
            # debug("{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
            # debug(filename)
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "papara_alignment.extended",
                         "-t", "place_resolve.tre",
                         "-p", "1",
                         "-n", "{}".format(self.date)])
        os.chdir(cwd)
        self._full_tree_est = 1  # Do we use these somewhere?

    def calculate_bootstrap(self):
        """calculate bootstrap and consensus trees
        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        """
        # debug("1")
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
        """This removes items from aln, and tree, if they were added to the blacklist.
        seq that were not added because they were similar to the one being removed here, are being lost
        """
        # that should not be a major issue though.
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            if gi_id in self.blacklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, gi is part of blacklist"
        self.data.reconcile()
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
                        # debug(filename.split("/")[-1])
                        # debug(self.blast_subdir)
                        if not os.path.exists("{}/previous_run".format(self.workdir)):
                            # debug('{}/previous_run/'.format(self.workdir))
                            os.makedirs('{}/previous_run/'.format(self.workdir))
                        # debug( "{}/{}".format(self.workdir, filename.split("/")[-1]))
                        # debug(filename)
                        os.rename(filename, "{}/{}".format(self.workdir, filename.split("/")[-1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file),
                          "{}/previous_run/newseqs.fasta".format(self.workdir))
                # try:
                self.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
                # except:
                #     self.data.write_labelled(label='^ot:ottTaxonName', gi_id=True)
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {}  # Wipe for next run
                self.new_seqs_otu_id = {}
                # self.newseqsgi = []  # Never replace it, is used to filter already added gis!!!
                self.repeat = 1
                # self.query_dict = {}  # clean up for next round
            else:
                if _VERBOSE:
                    sys.stdout.write("No new sequences after filtering.\n")
                self.repeat = 0
        else:
            # if glob.glob("/".join([self.workdir, "RAxML_bootstrap.all*"])):
            if _VERBOSE:
                sys.stdout.write("No new sequences found.\n")
            self.repeat = 0
            self.calculate_bootstrap()
            # else:
            #     self.repeat = 1
        self.reset_markers()
        self.data.dump()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))

    def write_unpl_lblastdb(self, path_to_local_seq):
        """Adds local sequences into a  local blast database, which then can be used to blast aln seq against it
        and adds sequences that were found to be similar to input.
        If this option is used, it queries against local database first and only in "2" round
        it goes back to blasting against GenBank"""
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
        # gi_counter = 1
        # add assert that tests that every file is a fasta file in the folder
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
            # that multiple seq fasta file can be read in as well
            for i in xrange(0, len(gi_id_l)):
                key = gi_id_l[i].replace(">", "")
                count = count + 1
                seq = seq_l[i]
                # debug(key)
                # debug(seq)
                write_blast_files(self.workdir, key, seq, db=True, fn="local_unpubl_seq")
        # general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
        # debug(cmd1)
        os.system(cmd1)


    def make_otu_dict_entry_unpubl(self, key):
        """adds the local unpublished data to the otu_dict.
        Information are retrieved from the additional json file/self.local_otu_json.
        I make up accession numbers....

        """
        debug("make_otu_dict_entry_unpubl")
        # debug(key)
        gi_counter = 1
        if key not in self.data.gi_dict.keys():
            # debug("key is new")
            # numbers starting with 0000 are unpublished data
            self.data.gi_dict[key] = {'accession': "000000{}".format(gi_counter), 'title': "unpublished", 'localID': key[7:]}
            gi_counter += 1
            # self.data.otu_dict[key] = {}
            # self.data.otu_dict[key]['^ncbi:gi'] = key
            # self.data.otu_dict[key]['^ncbi:accession'] = self.data.gi_dict[key]['accession']
            # self.data.otu_dict[key]['^user:TaxonName'] = self.data.gi_dict[key]['localID']
            # self.data.otu_dict[key]['^ncbi:title'] = self.data.gi_dict[key]['title']
            # local_id = self.data.gi_dict[key]['localID']
            # key2 = "otu{}".format(local_id)
            # self.data.otu_dict[key]['^ot:ottTaxonName'] = self.local_otu_json[key2]['^ot:ottTaxonName']
            # self.data.otu_dict[key]['^ncbi:taxon'] = self.local_otu_json[key2]['^ncbi:taxon']
            # self.data.otu_dict[key]['^ot:ottId'] = self.local_otu_json[key2]['^ot:ottId']
            # self.data.otu_dict[key]['^physcraper:status'] = "local seq"
            # self.data.otu_dict[key]['^physcraper:last_blasted'] = "1800/01/01"
            # self.ids.get_rank_info(taxon_name=self.data.otu_dict[key]['^ot:ottTaxonName'])
        else: 
            # debug("add new k,v - pairs")
            # debug(self.data.gi_dict[key])
            self.data.gi_dict[key].update([('accession', "000000{}".format(gi_counter)), ('title', "unpublished"), ('localID', key[7:])])


###############################

class FilterBlast(PhyscraperScrape):
    """Takes the Physcraper Superclass and filters the ncbi blast results to only include a subset of the sequences.

    They can be filtered by number or by rank and number. This can be useful for non-population-level studies,
    e.g. analyses which require a single representative per taxon (threshold = 1)
    or to check the monophyly of taxa without having to deal with over-representation of few taxa (e.g. threshold = 4,
    which allows to get a good overview of what is available without having some taxa being represented
    by high numbers of sequences). The second option (downtorank) allows to filter according to taxonomic levels,
    e.g. getting a number of representative sequences for a genus or lineage.
    This can also be used to not have to deal with subspecies.
    """

    def __init__(self, data_obj, ids_obj, settings=None):
        super(FilterBlast, self).__init__(data_obj, ids_obj)
        debug("start derived class init")
        # self.workdir = data_obj.workdir
        # self.logfile = "{}/logfile".format(self.workdir)
        # self.data = data_obj
        # self.ids = ids_obj
        # self.config = self.ids.config  # this is already part of .ids, or not? Information are doubled.
        # self.new_seqs = {}
        # self.new_seqs_otu_id = {}
        # self.otu_by_gi = {}  # What was this intended for?
        # self._to_be_pruned = []  # What is this used for?
        # self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        # self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        # self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        # if not os.path.exists(self.workdir):
        #     os.makedirs(self.workdir)
        # self.newseqs_file = "tmp.fasta"
        # self.date = str(datetime.date.today())
        # self.repeat = 1
        # self.reset_markers()
        # self.gi_list_mrca = []  # used for local blast to limit seq to seq of interest
        # if self.config.blast_loc == 'local' and len(self.gi_list_mrca) == 0:
        #     self.gi_list_mrca = self.get_all_gi_mrca()
        # self.unpublished = False
        # self.path_to_local_seq = False
        # self.local_otu_json = None
        # self.ncbi_parser = ncbi_data_parser.parser(names_file=self.config.ncbi_parser_names_fn, nodes_file=self.config.ncbi_parser_nodes_fn)

        # self.query_dict = {}  # for local blast txt files

        # additional things that are needed for the filtering process
        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}
        self.downtorank = None
        # self.localblast = False
        # self.local_otu_json = None
        # self.not_added = []
        if settings is not None:
            self.blacklist = settings.blacklist
        else:
            self.blacklist = []
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,",  "local"]

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
        # debug(self.data.otu_dict)
        self.sp_d = {}
        for key in self.data.otu_dict:
            # debug(key)
            # debug(key['^physcraper:status'])
            if self.data.otu_dict[key]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug(self.downtorank)
                spn = self.ids.find_name(sp_dict=self.data.otu_dict[key])
                if spn is None:
                    debug("spn is None")
                    # debug(self.data.otu_dict[key])
                    gi_id = self.data.otu_dict[key]['^ncbi:gi']
                    # debug(gi_id)
                    # debug(type(gi_id))
                    # spn = self.ids.get_rank_info(gi_id=gi_id)
                    spn = self.ids.find_name(gi=gi_id)
                    if spn is None:
                        debug("something is going wrong!Check species name")
                        sys.stderr.write("{} has no corresponding spn! Check what is wrong!".format(key))
                if self.downtorank is not None:
                    spn = str(spn).replace(" ", "_")
                    # tax_id = self.ids.get_ncbiid_from_tax_name(spn)
                    tax_id = self.ids.ncbi_parser.get_id_from_name(spn)
                    downtorank_id = self.ids.ncbi_parser.get_downtorank_id(tax_id, self.downtorank)
                    downtorank_name = self.ids.ncbi_parser.get_name_from_id(downtorank_id)

                    # if spn not in self.ids.otu_rank.keys():
                    #     self.ids.get_rank_info(taxon_name=spn)
                    #     spn = str(spn).replace(" ", "_")
                    # # if spn in self.ids.otu_rank.keys():
                    # lineage2ranks = self.ids.otu_rank[spn]["rank"]
                    # # else:
                    # #     self.ids.get_rank_info(taxon_name=spn)
                    # #     # debug(self.ids.otu_rank.keys())
                    # #     lineage2ranks = self.ids.otu_rank[str(spn).replace(" ", "_")]["rank"]
                    # #     # debug(lineage2ranks)s
                    # ncbi = NCBITaxa()
                    # for key_rank, val in lineage2ranks.iteritems():
                    #     if val == downtorank:
                    #         tax_id = key_rank
                    #         value_d = ncbi.get_taxid_translator([tax_id])
                    #         spn = value_d[int(tax_id)]
                    spn = downtorank_name
                if spn in self.sp_d:
                    self.sp_d[spn].append(self.data.otu_dict[key])
                else:
                    self.sp_d[spn] = [self.data.otu_dict[key]]
        return self.sp_d

    def make_sp_seq_dict(self):
        """Uses the sp_d to make a dict with species names as key1, key2 is gi/sp.name and value is seq: return sp_seq_d.

        This is used to select representatives during the filtering step, where it
        selects how many sequences per species to keep in the alignment. It will only contain sp that were not removed
        in an earlier cycle of the program.
        """
        # Note: has test, test_sp_seq_d.py, runs
        debug("make_sp_seq_dict")
        for key in self.sp_d:
            # loop to populate dict. key1 = sp name, key2= gi number, value = seq,
            # number of items in key2 will be filtered according to threshold and already present seq
            # debug(key)
            seq_d = {}
            for gi_id in self.sp_d[key]:
                # following if statement should not be necessary as it is already filtered in the step before.
                # I leave it in for now.
                if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    # I am using the next if to delimit which seq where already present from an earlier run,
                    # they will get a sp name (str), in order to distinguish them from newly found seq,
                    # which will have the gi (int). This differentiation is needed in the filtering blast step.
                    if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                        user_name = self.ids.find_name(sp_dict=gi_id)
                        for user_name_aln, seq in self.data.aln.items():
                            otu_dict_label = self.ids.find_name(sp_dict=self.data.otu_dict[user_name_aln.label])
                            if user_name == otu_dict_label:
                                seq = seq.symbols_as_string().replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[user_name_aln.label] = seq
                    else:
                        if "^ncbi:gi" in gi_id:  # this should not be needed: all new blast seq have gi
                            # debug("gi in gi_id")
                            # gi_num = int(gi_id['^ncbi:gi'])
                            gi_num = gi_id['^ncbi:gi']
                            if gi_num in self.new_seqs.keys():
                                seq = self.new_seqs[gi_num]
                                seq = seq.replace("-", "")
                                seq = seq.replace("?", "")
                                seq_d[gi_num] = seq
                    self.sp_seq_d[key] = seq_d
        # debug(self.sp_seq_d)
        return

    def select_seq_by_local_blast(self, seq_d, fn, treshold, count):
        """Selects number of sequences from local_blast to fill up to the threshold. It returns a filtered_seq dictionary.

        It will only include species which have a blast score of mean plus/minus sd.
        Is used after read_local_blast.
        """
        # Note: has test,test_select_seq_by_local_blast.py, runs
        debug("select_seq_by_local_blast")
        seq_blast_score = read_local_blast(self.workdir, seq_d, fn)
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
        # debug(self.filtered_seq)
        return self.filtered_seq

    def select_seq_by_length(self, taxon_id, treshold, count):
        """Select new sequences by length instead of by score values.
        """
        debug("select_seq_by_length")
        max_len = max(self.sp_seq_d[taxon_id].values())
        # !!! sometimes the only seq in seq_w_maxlen is the original seq,
        # then this is the one to be added, but it will be removed,
        # later as it is no new seq! thus no new seq for that species is added
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
        assert random_seq_ofsp is not None
        if random_seq_ofsp is not None:
            for key in random_seq_ofsp.keys():
                # debug(key)
                self.filtered_seq[key] = random_seq_ofsp[key]

    def add_all(self, key):
        """Add all seq to filtered_dict as the number of sequences is smaller than the threshold value.
        """
        # Note: has test, test_add_all.py: runs
        debug('add_all')
        for gi_id in self.sp_d[key]:
            if '^physcraper:status' in gi_id:
                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if gi_id['^physcraper:last_blasted'] == '1800/01/01':
                        gi_num = gi_id['^ncbi:gi']
                        # debug(gi_num)
                        seq = self.new_seqs[gi_num]
                        self.filtered_seq[gi_num] = seq
        return self.filtered_seq

    def get_name_for_blastfiles(self, key):
        """ Gets the name which is needed to write the blast files in 'loop_for_write_blast files'.

        The name needs to be retrieved before the actual loop starts. I use the taxonomic names here,
        as this is the measure of which information goes into which local filter blast database.
        The function is only used within 'loop_for_write_blast_files' to generate the filenames.
        """
        nametoreturn = None
        # loop needs to happen before the other one, as we need nametoreturn in second:
        for gi_id in self.sp_d[key]:
            # this if should not be necessary
            spn_name = self.ids.find_name(sp_dict=gi_id)
            if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                # if spn_name is None:
                #     # debug(key)
                #     spn_name = self.ids.get_rank_info(taxon_name=key)
                #     # gi_id['^ot:ottTaxonName'] = spn_name
                # spn_name = spn_name.replace(" ", "_")
                for spn_name_aln, seq in self.data.aln.items():
                    otu_dict_name = self.ids.find_name(sp_dict=self.data.otu_dict[spn_name_aln.label])
                    if spn_name == otu_dict_name:
                        nametoreturn = spn_name_aln.label
            # the next lines where added because the test was breaking,
            # needs thorough testing if it not breaks something else now.
            assert spn_name is not None  # assert instead of if
            if nametoreturn is None and spn_name is not None:
                nametoreturn = spn_name.replace(" ", "_")
            # debug(nametoreturn)
            # else:
            #     debug("do something?")
        return nametoreturn

    def loop_for_write_blast_files(self, key):
        """This loop is needed to be able to write the local blast files for the filtering step correctly.

        Function returns a filename for the filter blast, which were generated with 'get_name_for_blastfiles'.
        """
        # Note: has test,test_loop_for_blast.py: runs
        debug("length of sp_d key")
        # debug(len(self.sp_d[key]))
        nametoreturn = self.get_name_for_blastfiles(key)
        for gi_id in self.sp_d[key]:
            # debug("in writing file for-loop")
            # this if should not be necessary
            if '^physcraper:status' in gi_id and gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                # debug("generate files used for blast")
                if gi_id['^physcraper:last_blasted'] != '1800/01/01':  # old seq
                    spn_name = self.ids.find_name(sp_dict=gi_id)
                    for spn_name_aln, seq in self.data.aln.items():
                        otu_dict_name = self.ids.find_name(sp_dict=self.data.otu_dict[spn_name_aln.label])
                        if spn_name == otu_dict_name:
                            filename = nametoreturn
                            # filename = spn_name_aln.label
                            if self.downtorank is not None:
                                filename = key
                            # print(filename, seq)
                            write_blast_files(self.workdir, filename, seq)
                else:
                    # debug("make gilist as local blast database")
                    if "^ncbi:gi" in gi_id:
                        gi_num = int(gi_id['^ncbi:gi'])
                        # debug(gi_num)
                        # if selectby == "blast":  # should be obsolete now?!
                        # debug("new")
                        file_present = False
                        # debug(gi_num)
                        if gi_num in self.new_seqs.keys():
                            file_present = True
                        if file_present:  # short for if file_present == True
                            if '^physcraper:status' in gi_id:
                                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                                    filename = gi_num
                                    # debug("write seq to db")
                                    # debug(nametoreturn)
                                    seq = self.sp_seq_d[key][gi_num]
                                    if self.downtorank is not None:
                                        filename = key
                                        nametoreturn = key
                                    # debug(filename)
                                    write_blast_files(self.workdir, filename, seq, db=True, fn=nametoreturn)
                                    # blastfile_taxon_names[gi_num] = gi_num
                    namegi = key
        if self.downtorank is not None:
            nametoreturn = key
        if nametoreturn is None:
            nametoreturn = namegi
        assert nametoreturn is not None
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
        # this determines if a taxonomic  name / otu is already present in the aln and how many new seq. were found
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
        debug("length of sp_d")
        debug(len(self.sp_d))
        for taxon_id in self.sp_d:
            debug(taxon_id)
            count_dict = self.count_num_seq(taxon_id)
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            if len(self.sp_d[taxon_id]) <= treshold:  # add all stuff to self.filtered_seq[gi_n]
                self.add_all(taxon_id)
            elif len(self.sp_d[taxon_id]) > treshold:  # filter number of sequences
                debug("filter number of sequences")
                # debug(self.sp_seq_d[taxon_id].keys())
                if taxon_id in self.sp_seq_d.keys():
                    if selectby == "length":
                        # debug("{}, {}, {}".format(self.sp_seq_d[taxon_id], treshold, seq_present))
                        self.select_seq_by_length(taxon_id, treshold, seq_present)
#                        self.select_seq_by_length(self.sp_seq_d[taxon_id], treshold, seq_present)
                    elif selectby == "blast":
                        # if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
                        if 1 <= seq_present < treshold and count_dict["new_taxon"] is False and query_count != 0:
                            # debug("seq_present>0")
                            if query_count + seq_present > treshold:  # species is not new in alignment, make blast with existing seq
                                taxonfn = self.loop_for_write_blast_files(taxon_id)
                                # # next loop does not seem to be used
                                # for element in self.sp_d[taxon_id]:
                                #     if '^ot:ottTaxonName' in element:
                                #         blast_seq = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                #         blast_db = "{}".format(element['^ot:ottTaxonName']).replace(" ", "_")
                                if self.downtorank is not None:
                                    taxonfn = taxon_id
                                run_local_blast(self.workdir, taxonfn, taxonfn)
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], taxonfn, treshold, seq_present)
                            elif query_count + seq_present <= treshold:
                                self.add_all(taxon_id)
                        elif seq_present == 0 and count_dict["new_taxon"] is True and query_count >= 1:  # species is completely new in alignment
                            # debug("completely new taxon to blast")
                            # species is completely new in alignment, \
                            # make blast with random species
                            # debug(count_dict)
                            # debug(taxon_id)
                            # debug(self.sp_seq_d)
                            
                            # this causes to add some taxa twice to aln and phy!!! never use add_otu twice!
                            # for item in self.sp_d[taxon_id]:
                            #     if '^ncbi:gi' in item:
                            #         # if self.config.blast_loc == 'local':
                            #         #     localblast = True
                            #         # else:
                            #         #     localblast = False
                            #         self.data.add_otu(item['^ncbi:gi'], self.ids)
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
                            write_blast_files(self.workdir, str_db, seq)  # blast qguy
                            # debug("blast db new")
                            blast_db = self.sp_seq_d[taxon_id].keys()[1:]
                            # debug(blast_db)
                            for blast_key in blast_db:
                                seq = self.sp_seq_d[taxon_id][blast_key]
                                write_blast_files(self.workdir, blast_key, seq, db=True, fn=str_db)  # local db
                            # make local blast of sequences
                            run_local_blast(self.workdir, str_db, str_db)
                            if len(self.sp_seq_d[taxon_id]) + seq_present >= treshold:
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], str_db, treshold, seq_present)
                            elif len(self.sp_seq_d[taxon_id]) + seq_present < treshold:
                                self.add_all(taxon_id)
                else:
                    debug("taxon not in sp_seq_dict")
        return

    def replace_new_seq(self):
        """Replaces the self.new_seqs with the filtered_seq information.
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
        for gi_num in seq_not_added:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'not added, there are enough seq per sp in tre'
        for gi_num in keylist:
            for key in self.data.otu_dict.keys():
                # debug(self.data.otu_dict[key])
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        reduced_new_seqs_dic[key] = self.filtered_seq[gi_num]
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added, as representative of taxon'
                        # self.data.otu_dict[otu_id]['^ncbi:gi'] = gi_num
        # might not be necessary, I was missing some gi's
        # when I was replacing the original one.
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
        self.new_seqs = deepcopy(reduced_new_seqs)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic)  # TODO: key is not exactly same format as before in new_seqs_otu_id
        # the dict does not seem to be used for something...
        # set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return self.new_seqs

    def write_otu_info(self, downtorank=None):
        """Writes a table to file with taxon names and number of representatives - taxon_sampling.
        Write a file with all relevant GenBank info to file.
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


# ##################
def run_local_blast(workdir, blast_seq, blast_db, output=None):
    """Runs  a local blast to get measurement of differentiation between available sequences for the same taxon concept.

    The blast run will only be run if there are more sequences found than specified by the threshold value.
    When several sequences are found per taxon, blasts each seq against all other ones found for that taxon.
    The measure of differentiation will then be used to select a random representative from the taxon concept,
    but allows to exclude potential mis-identifications.
    In a later step (select_seq_by_local_blast) sequences will be selected based on the blast results generated here.
    """
    # Note: has test, runs -> test_run_local_blast.py
    debug("run_local_blast")
    debug(blast_seq)
    general_wd = os.getcwd()
    os.chdir(os.path.join(workdir, "blast"))
    out_fn = "{}_tobeblasted".format(str(blast_seq))
    cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(blast_seq)
    # debug("make local db")
    os.system(cmd1)
    if output is None:
        cmd2 = "blastn -query {} -db {}_db -out output_{}.xml -outfmt 5".format(out_fn, blast_db, out_fn)
    else:
        cmd2 = "blastn -query {} -db {}_db -out {} -outfmt 5".format(out_fn, blast_db, output)
    os.system(cmd2)
    # debug(cmd2)
    os.chdir(general_wd)


def calculate_mean_sd(hsp_scores):
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


def read_local_blast(workdir, seq_d, fn):
    """Reads the files of the local blast run and returns sequences below a value
    (within the std of the mean scores of the hsp.bit scores at the moment).

    (this is to make sure seqs chosen are representative of the taxon)
    """
    # Note: has test, runs: test_read_local_blast.py
    general_wd = os.getcwd()
    os.chdir(os.path.join(workdir, "blast"))
    output_blast = "output_{}_tobeblasted.xml".format(fn)
    xml_file = open(output_blast)
    os.chdir(general_wd)
    blast_out = NCBIXML.parse(xml_file)
    hsp_scores = {}

    tries = 5
    for i in range(tries):
        try:
            for record in blast_out:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        gi_id = alignment.title.split(" ")[1]
                        if gi_id.isdigit():
                            gi_id = int(gi_id)
                        # except:
                        #     gi_id = gi_id
                        # # debug(gi_id)
                        hsp_scores[gi_id] = {'hsp.bits': hsp.bits, 'hsp.score': hsp.score,
                                             'alignment.length': alignment.length, 'hsp.expect': hsp.expect}
        except ValueError:
            debug("rebuild the local blast db and try again")
            sys.stderr.write("{} blast file has a problem. Redo running it".format(fn))
            general_wd = os.getcwd()
            os.chdir(os.path.join(workdir, "blast"))
            # os.remove("{}_db.*".format(fn))
            subprocess.call(["rm", "{}_db.*".format(fn)])
            cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(fn)
            os.system(cmd1)
            cmd2 = "blastn -query {}_tobeblasted -db {}_db -out output_{}tobeblasted.xml -outfmt 5".format(fn, fn, fn)
            os.system(cmd2)
            os.chdir(general_wd)
            if i < tries - 1:  # i is zero indexed
                continue
            else:
                # debug("im going to raise")
                raise
        break
    # make values to select for blast search, calculate standard deviation,mean
    mean_sd = calculate_mean_sd(hsp_scores)
    # select which sequences to use
    seq_blast_score = {}
    for gi_id in hsp_scores:  # use only seq that are similar to mean plus minus sd
        # print(gi_id, hsp_scores[gi_id]['hsp.bits'])
        if (hsp_scores[gi_id]['hsp.bits'] >= mean_sd['mean'] - mean_sd['sd']) & \
                (hsp_scores[gi_id]['hsp.bits'] <= mean_sd['mean'] + mean_sd['sd']):
            if gi_id in seq_d:
                seq_blast_score[gi_id] = seq_d[gi_id]
    return seq_blast_score


def write_blast_files(workdir, file_name, seq, db=False, fn=None):
    """Writes local blast files which will be read by run_local_blast.
    """
    debug("writing files")
    # debug(file_name)
    if not os.path.exists("{}/blast".format(workdir)):
        os.makedirs("{}/blast/".format(workdir))
    if db:
        fnw = "{}/blast/{}_db".format(workdir, fn)
        fi_o = open(fnw, 'a')
    else:
        fnw = "{}/blast/{}_tobeblasted".format(workdir, file_name)
        fi_o = open(fnw, 'w')
    # debug(fnw)
    fi_o.write(">{}\n".format(file_name))
    fi_o.write("{}\n".format(str(seq).replace("-", "")))
    fi_o.close()


def get_ncbi_tax_id(handle):
    """Get the taxon ID from ncbi.
    """
    gb_list = handle[0]['GBSeq_feature-table'][0]['GBFeature_quals']
    for item in gb_list:
        if item[u'GBQualifier_name'] == 'db_xref':
            # debug(item[u'GBQualifier_value'])
            if item[u'GBQualifier_value'][:5] == 'taxon':
                ncbi_taxonid = int(item[u'GBQualifier_value'][6:])
                break
            else:
                continue
    return ncbi_taxonid


def get_ncbi_tax_name(handle):
    """Get the sp name from ncbi.
    """
    gb_list = handle[0]['GBSeq_feature-table'][0]['GBFeature_quals']
    for item in gb_list:
        if item[u'GBQualifier_name'] == 'organism':
            ncbi_sp = str(item[u'GBQualifier_value'])
    return ncbi_sp
