#!/usr/bin/env python
"""Physcraper module"""
import sys
import re
import os
import subprocess
import time
import datetime
import glob
import json
import unicodedata
import configparser
import pickle
import collections
import inspect
import random
import logging
import numpy
import collections 
import inspect
from copy import deepcopy
from ete2 import NCBITaxa
from urllib2 import URLError
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from dendropy import Tree,\
                  DnaCharacterMatrix,\
                  DataSet,\
                  datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree,\
                              extract_tree_nexson,\
                              get_subtree_otus,\
                              extract_otu_nexson,\
                              PhyloSchema
from peyotl.api import APIWrapper
import physcraper.AWSWWW as AWSWWW


_DEBUG = 1
_DEBUG_MK = 1

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


class ConfigObj(object):
    """Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store."""
    def __init__(self, configfi):
        if _DEBUG:
            sys.stderr.write("Building config object\n")
        debug(configfi)
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
        if not os.path.isfile(self.ncbi_dmp):
            os.system("rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz {}.gz".format(self.config.ncbi_dmp))
            os.system("tar -xzvf taxonomy/gi_taxid_nucl.dmp.gz")
            self.ncbi_dmp = "taxonomy/gi_taxid_nucl.dmp.gz"
        self.phylesystem_loc = config['phylesystem']['location']
        assert(self.phylesystem_loc in ['local', 'api'])
        self.ott_ncbi = config['taxonomy']['ott_ncbi']
        assert os.path.isfile(self.ott_ncbi)
        self.id_pickle = os.path.abspath(config['taxonomy']['id_pickle'])#TODO what is theis doing again?
        self.email = config['blast']['Entrez.email']
        assert '@' in self.email
        self.blast_loc = config['blast']['location']
        self.num_threads = config['blast'].get('num_threads')
        assert self.blast_loc in ['local', 'remote']
        if self.blast_loc =='local':
            self.blastdb = config['blast']['localblastdb']
            self.url_base = None
        if self.blast_loc =='remote':
            self.url_base = config['blast'].get('url_base')
        if _DEBUG:
            sys.stderr.write("{}\n".format(self.email))
            if self.blast_loc =='remote':
                sys.stderr.write("url base = {}\n".format(self.url_base))
            sys.stderr.write("{}\n".format(self.blast_loc))
            if self.blast_loc =='local':
                sys.stderr.write("local blast db {}\n".format(self.blastdb))


#ATT is a dumb acronym for Alignment Tree Taxa object
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
#        url = "http://treebase.org/treebase-web/phylows/study/TB2:S{}?format=nexml".format(tb_id)
        if _DEBUG:
             sys.stderr.write(url +"\n")
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
    #TODO CHECK ARGS
    assert isinstance(aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH
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
    newick = newick.replace(" ", "_") #UGH Very heavy handed, need to make sure happens on alignement side as well.
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
            sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignement. This may indicate a mismatch between tree and alignement\n".format(tax.label))
   #need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    workdir = os.path.abspath(workdir)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)
    #newick should be bare, but alignement should be DNACharacterMatrix


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
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH
    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    assert tre.taxon_namespace is aln.taxon_namespace
    otu_newick = tre.as_string(schema=schema_trf)
    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        otu_dict = json.load(open(otu_json,"r"))
        ott_ids = [otu_dict[otu].get(u'^ot:ottId',) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
        debug(ott_mrca)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)


#def otu_tiplabel(aln):
#    """makes a dictionary with tiplabels of aln as key and
#    removes 'n' from tiplabels that start with a number as value.
#    """
#    otu_taxonlabel_problem = {}
#    otu_taxonlabel_problem[key.label] = newname
#    return otu_taxonlabel_problem


def standardize_label(item):
        """try to make names unicode
        """
        # Note: has test, runs -> test_edit_dict_key.py
        debug("standardize label")
        item_edit = item.replace("-", "_")
        item_edit = item_edit.replace("'", "_")
        item_edit = item_edit.replace("/", "_")
        return item_edit

def OtuJsonDict(id_to_spn, id_dict):
    """Make otu json dict, which is also produced within the openTreeLife-query
    reads input file into the var spInfo, tranaltes using an IdDict object
    using web for NCBI"""

    with open(id_to_spn, mode='r') as idtospn:
        reader = csv.reader(idtospn)
        spInfo = dict((rows[0], rows[1]) for rows in reader)
    ncbi = NCBITaxa()    
    spInfoDict = {}
    for item in spInfo:
        spn = spInfo[item].replace("_", " ")
        name2taxid = ncbi.get_name_translator([spn])
        clean_lab = standardize_label(item)
        otuid = "otu{}".format(clean_lab)
        if len(name2taxid.items())>=1:
            ncbiid = name2taxid.items()[0][1][0]
            ott = id_dict.ncbi_to_ott[ncbiid]
            spn = id_dict.ott_to_name[ott]
            get_info = {'^ncbiID': ncbiid, '^ot:ottTaxonName': spn, '^ot:ottId': ott, '^ot:originalLabel': item, '^user:TaxonName': spInfo[item],  '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01" }  
            spInfoDict[otuid] = get_info
        else:
            
            spInfoDict[otuid] = {'^user:TaxonName': spInfo[item], '^ot:originalLabel': item, '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01"}
    return  spInfoDict 

class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.
    Hypothetically, all teh keys in the  otu_dict shopuld be cealn
    """
    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None, taxon_namespace=None):
        #TODO add assertions that inputs are correct type!!!
        self.aln = alignment
        if schema == None:
            self.tre = Tree.get(data=newick,
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        else:
            self.tre = Tree.get(data=newick,
                                schema=schema,
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        assert(self.tre.taxon_namespace is self.aln.taxon_namespace)
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
        self.tre = Tree.get(data=newick,
                            schema="newick",
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        assert isinstance(otu_dict, dict)
        self.otu_dict = otu_dict
        self.ps_otu = 1 #iterator for new otu IDs
#        self.otu_taxonlabel_problem = otu_tiplabel(alignment) #Check back - original labels are not stored?
        self._reconcile_names()
        self.workdir = os.path.abspath(workdir) #TODO - is this where the workdir should live?
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca)
        self.ott_mrca = ingroup_mrca
        self.orig_seqlen = [] #FIXME
        self.gi_dict = {}
        self.orig_aln = alignment
        self.orig_newick = newick

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
            try:
                self.otu_dict[tax.label]['^physcraper:status'] = "deleted in name reconciliation"
            except:
                sys.stderr.write("lost taxon {} in name reconcilliation".format(tax.label))
           #     self.otu_dict[self.otu_taxonlabel_problem[tax.label]]['^physcraper:status'] = "deleted in name reconciliation" #ejm look
            self.aln.taxon_namespace.remove_taxon(tax)
        assert(self.aln.taxon_namespace==self.tre.taxon_namespace)
        reverse_otu_dict = {}
        for tax in aln.taxon_namespace:
            found_label = 0
            match = re.match("'n[0-9]{1,3}", tax.label)
            newname = ''
            if match:
                newname = tax.label[2:]
                newname = newname[:-1]
            for otu in self.otu_dict:
                if self.otu_dict[otu]['^ot:originalLabel'] == tax.label or self.otu_dict[otu]['^ot:originalLabel'] == newname:
                    tax.label == otu
                    found_label = 1
            if found_label == 0:
                sys.stderr.write("could not match tiplabel {} or {} to an OTU".format(tax.label, newname))
        #TODO - make sure all taxon labels are unique OTU ids.


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
            debug(prune)
            self.aln.remove_sequences(prune)
            self.tre.prune_taxa(prune)
            self.tre.prune_taxa_with_labels(prune)#sometimes it does not delete it with the statement before. Tried to figure out why, have no clue yet.
            #self.aln.taxon_namespace.remove_taxon_label(tax)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short step due to sequence shorter than {}\n".format(min_seqlen))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in prune short"
            self.aln.taxon_namespace.remove_taxon_label(tax.label) # raises error if not found, instead of remove_taxon
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        self.reconcile()


    def reconcile(self, seq_len_perc=0.75): #
        """all missing data seqs are sneaking in, but from where?!"""
        # not only used in the beginning...is used to remove sequences that are shorter than 75%
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        debug("reconcile")
        prune = []
        debug(self.orig_seqlen)
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        debug(avg_seqlen)
        seq_len_cutoff = avg_seqlen*seq_len_perc
        debug(seq_len_cutoff)
        for tax, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(tax)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in reconcilation step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict.get(tax.label).get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            debug(tax)
            debug(tax.label)
            debug(self.otu_dict[tax.label])

            self.otu_dict[tax.label]['^physcraper:status'] = "deleted in reconcile" # unnecessary? gehts overwriten in next line
            self.remove_taxa_aln_tre(tax.label)
        aln_ids = set()
        self.edit_dict_key()
        for tax in self.aln:
            aln_ids.add(tax.label)
        assert aln_ids.issubset(self.otu_dict.keys())
        treed_taxa = set()
        orphaned_leafs = set()
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        
        ## here leaf_nodes have taxa that were dropped before. Why do we have this anyways?  
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
            if leaf.taxon.label not in aln_ids:
                self.otu_dict[leaf.taxon.label]['^physcraper:status'] = "deleted due to presence in tree but not aln. ?!"
                orphaned_leafs.add(leaf)
                #TODO fiure out why sometimes one of them works ano not the other and vice versa
                self.tre.prune_taxa([leaf])
                self.tre.prune_taxa_with_labels([leaf.taxon.label])
                self.tre.prune_taxa_with_labels([leaf])
                treed_taxa.remove(leaf.taxon.label)
                #debug(self.otu_taxonlabel_problem.keys())
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



    def trim(self, taxon_missingness = 0.75):
        '''cuts off ends of alignment, maintaining similar to original seq len
        Important bc other while whole chromosomes get dragged in!'''
        debug('in trim')
        seqlen = len(self.aln[0])
        for tax in self.aln:
            if len(self.aln[tax]) != seqlen:
                sys.stderr.write("can't trim un-aligned inputs, moving on")
                return
        start = 0
        stop = seqlen
        cutoff = len(self.aln) *  taxon_missingness
        for i in range(seqlen):
            counts = {'?':0, '-':0}
            debug(counts)
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?']+counts['-'] <= cutoff: #first ok column
                start = i
                break
        for i in range(seqlen-1, 0, -1):
            counts = {'?':0, '-':0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?']+counts['-'] <= cutoff:
                stop = i
                break
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        sys.stdout.write("trimmed alignement ends to < {} missing taxa, start {}, stop {}\n".format(taxon_missingness, start, stop))
        return


    def add_otu(self, gi, ids_obj):
        """generates an otu_id for new sequences and adds them into the otu_dict.
        Needs to be passed an IdDict to do the mapping"""
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1

        ncbi_id = ids_obj.map_gi_ncbi(gi) #check that try/execpt is working here
        #TODO do we need rank info here?
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi]['title']
        self.otu_dict[otu_id]['^ncbi:taxon'] = ncbi_id
        self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ncbi_id)
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = ids_obj.ott_to_name.get(self.otu_dict[otu_id]['^ot:ottId'])
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01" #1800 = never blasted; 1900 = blasted 1x, not added, this century = blasted and added 
        if _DEBUG >= 2:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi, otu_id))
        return otu_id



    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """Papara is finicky about trees and needs phylip, this writes out needed files for papara
        (except query sequences)"""
        #CAN I even evaulte things in the function definitions?
        debug('write papara files')
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")#Papara is diffffffficult about root
        tmptre = tmptre.replace("'", "_")
        fi = open("{}/{}".format(self.workdir, treefilename), "w")
        fi.write(tmptre)
        fi.close()
        self.aln.write(path="{}/{}".format(self.workdir, alnfilename), schema="phylip")
        # hack to remove illegal characters in files
        self.rewrite_files(inputfn="random_resolve.tre")
        self.rewrite_files(inputfn="aln_ott.phy")


    def write_files(self, treepath="physcraper.tre", treeschema="newick", alnpath="physcraper.fas", alnschema="fasta"):
        """Outputs both the streaming files and a ditechecked"""
        #First write rich annotation json file with everything needed for later?
        debug("write_files")
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema, unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)
        # hack to remove illegal characters in files

    def write_labelled(self, label, treepath="labelled.tre", alnpath="labelled.fas"):
        """output tree and alignement with human readable labels
        Jumps through abunch of hoops to make labels unique.
        NOT MEMORY EFFICIENT AT ALL"""
        debug("write labelled files")
        assert label in ['^ot:ottTaxonName', 'user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        tmp_fasta = self.aln.as_string(schema="fasta")
        tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta,
                                         schema="fasta")
        new_names = set()  
        for taxon in tmp_tre.taxon_namespace:
            debug(taxon)
            ### here the double names of the labelled tre files are generated.
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label == None:
                if otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_".format(otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(otu_dict[taxon.label].get("^ncbi:taxon", "unk"), otu_dict[taxon.label].get('^ot:ottTaxonName', "unk"))
            new_label = new_label.replace(' ', '_')
            if new_label in new_names:
                new_label = " ".join([new_label, taxon.label])
                new_names.add(new_label)
            taxon.label = new_label
        tmp_tre.write(path="{}/{}".format(self.workdir, treepath),
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path="{}/{}".format(self.workdir, alnpath),
                      schema="fasta")



    def write_otus(self, filename, schema='table'):
        """Writes out OTU dict as json"""
        assert schema in ['table', 'json']
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)

    def remove_taxa_aln_tre(self, taxon_label):
        """removes taxa from aln, tre and otu_dict,
        takes a single taxon_label as input.
        """
        # note: has test, test_remove_taxa_aln_tre.py, runs, passes
        debug('remove_taxa_aln_tre')
        debug(taxon_label)
        debug(type(taxon_label))
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        debug(tax)
        debug(tax2)
        #assert tax==tax2
        debug(len(self.tre.taxon_namespace))
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label) # raises an error if label not found
            #the first prune does not remove it sometimes...
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            # self.tre.taxon_namespace.remove_taxon_label(taxon_label) # cannot happen here, as then the seq_dict_build_ crashes, if taxon was just added?
            debug(len(self.tre.taxon_namespace))
            self.otu_dict[tax.label]['^physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['^physcraper:status'] = "deleted, but it wasn't in the alignment..."


    def dump(self, filename = None ):
        if filename:
            ofi = open(filename, "wb" )
        else:
            ofi = open("{}/att_checkpoint.p".format(self.workdir), "wb" )
        pickle.dump(self, ofi)

#####################################
def get_nexson(study_id, phylesystem_loc):
    """Grabs nexson from phylesystem"""
    phy = PhylesystemAPI(get_from=phylesystem_loc)
    nexson = phy.get_study(study_id)['data']
    return  nexson

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
    mrca_node = tree_of_life.mrca(ott_ids=synth_tree_ott_ids, wrap_response=False)# need to fix wrap eventually
    try:
        tax_id = mrca_node[u'mrca'][u'taxon'][u'ott_id']
        sys.stdout.write('MRCA of sampled taxa is {}\n'.format(mrca_node[u'mrca'][u'taxon'][u'name']))
    except: #Hackaround for V2 apis
        tax_id = mrca_node[u'nearest_taxon_mrca_ott_id']
        sys.stdout.write('MRCA of sampled taxa is {}\n'.format(mrca_node[ u'nearest_taxon_mrca_unique_name']))
    return tax_id


def get_ott_ids_from_otu_dict(otu_dict): #TODO put into data obj?
    """Get the ott ids from an otu dict object"""
    ott_ids = []
    for otu in otu_dict:
        try:
            ott_ids.append(otu['^ot:ottId'])
        except KeyError:
            pass
#####################################

class IdDicts(object):
    """Wraps up the annoying conversions"""#TODO - could - should be shared acrosss runs?! .... nooo.
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
        fi = open(config_obj.ott_ncbi) #TODO need to keep updated
        for lin in fi: #TODO This is insanely memory inefficient
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()
            assert len(self.ott_to_ncbi) > 0
            assert len(self.ncbi_to_ott) > 0
            assert len(self.ott_to_name) > 0
        fi.close()
        if os.path.isfile("{}/id_map.txt".format(workdir)): #todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]


    def get_rank_info(self, gi_id=False, taxon_name=False):
        """collect rank and linegae information from ncbi,
        used to delimit the sequences from blast,
        when you have a local blast database
        """
        ##I do this, to limit the results to the taxon group, as usually done directly during the blast search.
        # do that now with get_all_gi_mrca, this function is stil used somewhere else
        # debug("get_rank_info")
        Entrez.email = self.config.email
        if gi_id != False:
            debug("gi_id to tax_name")
            tries = 5
            for i in range(tries):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=gi_id, retmode="xml")
                except:
                    print(i)
                    if i < tries - 1: # i is zero indexed
                        continue
                    else:
                        raise
                break
            read_handle = Entrez.read(handle)[0]
            tax_name = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value']
            # debug(tax_name)
        else:
            tax_name = taxon_name.replace("_", " ")
        # debug(tax_name)
        if tax_name not in self.otu_rank.keys():
            debug("tax_name to rank")
            try:
                debug("try")
                tax_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
            except:
                debug("except")
                ncbi = NCBITaxa()
                tax_info = ncbi.get_name_translator([tax_name])
                debug(tax_info)
                if tax_info == {}:
                    print("Taxon name does not match any species name in ncbi. Check that the name is written correctly!")
                tax_id = tax_info.items()[0][1][0]
            ncbi = NCBITaxa()
            lineage = ncbi.get_lineage(tax_id)
            lineage2ranks = ncbi.get_rank(lineage)
            tax_name = str(tax_name).replace(" ", "_")
            self.otu_rank[tax_name] = {"taxon id": tax_id, "lineage": lineage, "rank": lineage2ranks}
        return tax_name


    def map_gi_ncbi(self, gi, get_rank = True):
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

    def dump(self, filename = None):
        if filename:
            ofi = open(filename, "wb" ) 
        else:
            ofi = open("{}/id_pickle.p".format(self.workdir,filename), "wb" ) 
        pickle.dump(self, ofi)
        


class PhyscraperScrape(object): #TODO do I wantto be able to instantiate this in a different way?!
    #set up needed variables as nones here?!
    #TODO better enforce ordering
    """This is the class that does the perpetual updating"""
    def __init__(self, data_obj, ids_obj):
        #todo check input types assert()
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
        self.date = str(datetime.date.today()) #Date of the run - may lag behind real date!
        self.repeat = 1
        self.reset_markers()

 #TODO is this the right place for this?
    def reset_markers(self):
        self._blasted = 0
        self._blast_read = 0
        self._identical_removed = 0
        self._query_seqs_written = 0
        self._query_seqs_aligned = 0
        self._query_seqs_placed = 0
        self._reconciled = 0
        self._full_tree_est = 0

    def run_blast(self): #TODO Should this be happening elsewhere?
        """generates the blast queries and saves them to xml"""
        debug("run_blast")
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            #TODO temp until I fix delete
            if otu_id in self.data.otu_dict:
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                if abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days) > 14: #TODO make configurable
                    equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi,
                                                                   last_blast,
                                                                   today)
                    query = seq.symbols_as_string().replace("-", "").replace("?", "")
                    xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
                    if not os.path.isfile(xml_fi):
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
                            #debug(blastcmd)
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
                                # debug(result_handle.read())
                            save_file = open(xml_fi, "w")
                            save_file.write(result_handle.read())
                            result_handle.close()
                            save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                       # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                       #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
                    else:
                        # changes date of blasted accordingly, if file is already present in the folder
                        if _DEBUG_MK == 1:
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
        self._blasted = 1


    def get_all_gi_mrca(self):
        """get all available gi numbers from Genbank for mrca.
        The list will be used to filter out sequences from the local Blast search,
        that do not belong to ingroup."""
        debug("get_all_gi_mrca")
        Entrez.email = self.config.email
        handle = Entrez.esearch(db="nucleotide", term="txid{}[Orgn]".format(self.mrca_ncbi), RetMax=1000000)
        records = Entrez.read(handle)
        id_list = records['IdList']
        id_list = [int(x) for x in id_list]
        return id_list


    def read_blast(self, blast_dir=None):
        """reads in and prcesses the blast xml files"""
        debug("read blast")
        if blast_dir:
            sys.stdout.write("blast dir {}".format(blast_dir))
            self.blast_subdir = os.path.abspath(blast_dir)
        sys.stdout.write("blast dir {}".format(self.blast_subdir))
        assert os.path.exists(self.blast_subdir)
        if not self._blasted:
            self.run_blast()
        if self.config.blast_loc == 'local': #because local db doens't have taxon info, needed to limit to group of interest
            gi_list_mrca = self.get_all_gi_mrca()
        for taxon in self.data.aln:
            debug("add blast seq to new seqs")
            xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    sys.stdout.write(".")
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    gi_id = int(alignment.title.split('|')[1])
                                    if gi_list_mrca and (gi not in gi_list_mrca):
                                        pass
                                    else:
                                       if gi_id not in self.data.gi_dict: #skip ones we already have (does it matter if these were deleted? No...)
                                                self.new_seqs[gi_id] = hsp.sbjct
                                                self.data.gi_dict[gi_id] = alignment.__dict__
                except ValueError:
                    sys.stderr.write("Problem reading {}, skipping\n".format(xml_fi))
        self.date = str(datetime.date.today())
        self._blast_read = 1

    # TODO this should go back in the class and should prune the tree


    def seq_dict_build(self, seq, label, seq_dict): #Sequence needs to be passed in as string.
        """takes a sequence, a label (the otu_id) and a dictionary and adds the
        sequence to the dict only if it is not a subsequence of a
        sequence already in the dict.
        If the new sequence is a super suquence of one in the dict, it
        removes that sequence and replaces it"""
        #TODO unify spp name somehow?
        if '^user:TaxonName' in self.data.otu_dict[label].keys():
            spn_of_label = self.data.otu_dict[label]['^user:TaxonName']
        elif '^ot:ottTaxonName' in self.data.otu_dict[label].keys():
            spn_of_label = self.data.otu_dict[label]['^ot:ottTaxonName']
        if spn_of_label is None:
            try:
                spn_of_label = self.data.otu_dict[label]['^ncbi:taxon']
            except:
                spn_of_label = None
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        for tax_lab in tax_list:
            if '^user:TaxonName' in self.data.otu_dict[tax_lab].keys():
                existing_taxa = self.data.otu_dict[tax_lab]['^user:TaxonName']
            elif '^ot:ottTaxonName' in self.data.otu_dict[tax_lab].keys():
                existing_taxa = self.data.otu_dict[tax_lab]['^ot:ottTaxonName']
                if existing_taxa is None:
                    existing_taxa = self.data.otu_dict[tax_lab]['^ncbi:taxon']
            existing_taxa = existing_taxa.replace("_", " ")
            spn_of_label = spn_of_label.replace("_", " ")
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(inc_seq) >= len(new_seq):
                if inc_seq.find(new_seq) != -1:
                    ## changed the code, that seq which are identical but belong to different species concepts, will be added
                    if spn_of_label != existing_taxa:
                        sys.stdout.write("seq {} is subsequence of {}, but different species concept\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; subsequence, but different species"
                        seq_dict[label] = seq
                        if _DEBUG_MK == 1:
                            print(spn_of_label, "and", existing_taxa, "subsequences, but different sp. concept")
                    else:
                        sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        if _DEBUG_MK == 1:
                            print(spn_of_label, " not added, subseq of ", existing_taxa)
                    return
            else:
                if new_seq.find(inc_seq) != -1:#this statement adds seq that are the same
                    if self.data.otu_dict[tax_lab].get('^physcraper:status') == "original":
                        # print("delete because its a superseq")

                        sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax_lab))

                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        if _DEBUG_MK == 1:
                            print(spn_of_label, " and ", existing_taxa, "added")
                        return
                    else:
                        # print("delete because its a superseq")
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added in place of {}".format(tax_lab)
                        if _DEBUG_MK == 1:
                            print(spn_of_label, "added, instead of ", existing_taxa)
                        return
        sys.stdout.write(".")
        if i%50 == 0:
            sys.stdout.write("\n")
        seq_dict[label] = seq
        return


    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        Does not test if they are identical to ones in the original alignment...."""
        debug("remove identical seqs")
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        #Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.data.orig_seqlen)/len(self.data.orig_seqlen) #HMMMMMMMM
        assert self.config.seq_len_perc <= 1
        seq_len_cutoff = avg_seqlen*self.config.seq_len_perc
        for gi, seq in self.new_seqs.items():
            if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                otu_id = self.data.add_otu(gi, self.ids)                   
                self.seq_dict_build(seq, otu_id, tmp_dict)
        for tax in old_seqs:
            try:
                del tmp_dict[tax]
            except KeyError:
                pass
        self.new_seqs_otu_id = tmp_dict # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from genbank, of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()


    def dump(self, filename = None):
        if filename:
            ofi = open(filename, "wb" )
        else:
            ofi = open("{}/scrape_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self,  ofi)
        #TODO... write as proper nexml?!  
         
    def write_query_seqs(self):
        """writes out the query sequence file"""
        debug("write query seq")
        if not self._blast_read:
            self.read_blast()
        self.newseqs_file = "{}.fasta".format(self.date)
        fi = open("{}/{}".format(self.workdir, self.newseqs_file), 'w')
        sys.stdout.write("writing out sequences\n")
        for otu_id in self.new_seqs_otu_id.keys():
            if otu_id not in self.data.aln: #new seqs only
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1
 
    def align_query_seqs(self, papara_runname="extended"):
        """runs papara on the tree, the alinment and the new query sequences"""
        cwd = os.getcwd()
        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        sys.stdout.write("aligning query sequences \n")
        ## hack around stupid characters for phylogen. tools
        self.data.edit_dict_key()
        ### note: sometimes there are still sp in any of the aln/tre and I still have not found out why sometimes the label is needed
        for tax_lab in self.data.aln.taxon_namespace:
            if tax_lab not in self.data.tre.taxon_namespace:
                debug("tax not in tre ")
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_lab in self.data.tre.taxon_namespace:
            if tax_lab not in self.data.aln.taxon_namespace:
                debug("tax not in aln ")
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_aln in self.data.aln.taxon_namespace:
            debug(tax_aln)
            if tax_aln not in self.data.tre.taxon_namespace:
                debug("tax not in tre ")
                self.data.remove_taxa_aln_tre(tax_aln.label)
        for tax_tre in self.data.tre.taxon_namespace:
            if tax_tre not in self.data.aln.taxon_namespace:
                debug("tax not in aln ")
                self.data.remove_taxa_aln_tre(tax_tre.label)
        ## hack for the alien taxa thing
        debug("added another reconcile step, maybe that helps with the S.glaber/alien taxa problem")
        self.data.reconcile()
        self.data.write_papara_files()
        os.chdir(self.workdir)#Clean up dir moving
        try:
            debug("I call papara")
            assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
            subprocess.call(["papara",
                             "-t", "random_resolve.tre",
                             "-s", "aln_ott.phy",
                             "-q", self.newseqs_file,
                             "-n", papara_runname]) #FIx directory ugliness
            sys.stdout.write("Papara done")
            self.data.rewrite_files(inputfn="papara_alignment.extended")
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
        debug(self.data.aln.taxon_namespace)
        self.data.aln.taxon_namespace.is_mutable = True #Was too strict...
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
        sys.stdout.write("placing query sequences \n")
        cwd = (os.getcwd())
        os.chdir(self.workdir)
        try:
            p1 = subprocess.call(["raxmlHPC", "-m", "GTRCAT",
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
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                              "-s", "papara_alignment.extended",
                              "-t", "place_resolve.tre",
                              "-p", "1",
                              "-n", "{}".format(self.date)])
        os.chdir(cwd)
        self._full_tree_est = 1


    def calculate_bootstrap(self):
        """calculate bootstrap and consensus trees
        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bottstrap stopping criteria
        """
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        # run bootstrap
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "aln_ott.phy",
                         "-t", "place_resolve.tre",
                         "-p", "1", "-b", "1", "-#", "autoMRE",
                         "-n", "bootstrap_{}".format(self.date)])
        # make bipartition tree
        # is the -f b command
        # -z specifies file with multiple trees
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "aln_ott.phy",
                         "-t", "RAxML_bestTree.{}".format(self.date),
                         "-p", "1", "-f", "b", "-z", "RAxML_bootstrap_{}".format(self.date),
                         "-n", "bipart_{}".format(self.date)])
        #strict consensus:
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-J", "STRICT",
                         "-z", "RAxML_bootstrap_{}".format(self.date),
                         "-n", "StrictCon{}".format(self.date)])
        #majority rule:
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-J", "MR",
                         "-z", "RAxML_bootstrap_{}".format(self.date),
                         "-n", "MR_{}".format(self.date)])
        #extended majority rule:
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-J", "MRE",
                         "-z", "RAxML_bootstrap_{}".format(self.date),
                         "-n", "EMR{}".format(self.date)])
        # rapid bootstrapping
        # -f a: is the command to do that
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-f", "a",
                         "-p", "1", "-x", "1", "-#", "autoMRE", "-s", "aln_ott.phy",
                         "-n", "rapidBS_{}".format(self.date)])


    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignment with the expanded ones"""
        ## first if should not be necessary, as this is in the wrapper
        if len(self.new_seqs) > 0:
            self.remove_identical_seqs()
            self.data.write_files() #should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:#TODO rename to something more intutitive
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
                os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
                if os.path.exists("{}/last_completed_update".format(self.workdir)):
                    os.rename(self.tmpfi, "{}/last_completed_update".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file), "{}/previous_run/newseqs.fasta".format(self.workdir))
                try:
                    self.data.write_labelled(label=  '^ot:ottTaxonName')
                except:
                    self.data.write_labelled(label='user:TaxonName')
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {} #Wipe for next run
                self.new_seqs_otu_id = {}
                self.repeat = 1
            else:
                sys.stdout.write("No new sequences after filtering.\n")
                self.repeat = 0
        else:
            sys.stdout.write("No new sequences found.\n")
            self.repeat = 0
            self.calculate_bootstrap()
        self.reset_markers()
        self.data.dump()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))


###############################

class FilterBlast(PhyscraperScrape):
    """takes the Physcraper Superclass and filters the ncbi blast results."""
    def __init__(self, data_obj, ids_obj):
        super(FilterBlast, self).__init__(data_obj, ids_obj)
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
        self.date = str(datetime.date.today())
        self.repeat = 1
        self.reset_markers()

        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}
        self.not_added = []
        self.seq_filter = ['deleted', 'subsequence,', 'not']

    def sp_dict(self, downtorank=None):
        """makes dict with species name as key and the corresponding seq
        information from aln and blast seq
        """
        # Note: has test, test_sp_d.py, runs
        self.downtorank = downtorank
        debug("make sp_dict")
        self.sp_d = {}
        for key in self.data.otu_dict:
            if downtorank != None:
                debug("downto is not None")
                ncbi = NCBITaxa()
                if '^ot:ottTaxonName' in self.data.otu_dict[key]:
                    tax_name = self.data.otu_dict[key]['^ot:ottTaxonName']
                    # ncbiID = ncbi.get_name_translator(tax_name)
                    ## does not always seem to work
                elif '^user:TaxonName' in self.data.otu_dict[key]:
                    tax_name = self.data.otu_dict[key]['^user:TaxonName']
                    # ncbiID = ncbi.get_name_translator(tax_name)
                # elif '^ncbi:taxon' in self.data.otu_dict[key]:
                #     ncbiID = self.data.otu_dict[key]['^ncbi:taxon']
                # elif '^ncbiID' in self.data.otu_dict[key]:
                #     ncbiID = self.data.otu_dict[key]['^ncbiID']
                else:
                    debug("no taxon name provided! It will fail")
                debug(tax_name)
                self.ids.get_rank_info(taxon_name=tax_name)
                debug(self.ids.otu_rank.keys())
                lineage2ranks = self.ids.otu_rank[tax_name.replace(" ", "_")]["rank"]
                for key_rank, val in lineage2ranks.iteritems():
                    if val == downtorank:
                        tax_id = key_rank
                        value_d = ncbi.get_taxid_translator([tax_id])
                        value = value_d[int(tax_id)]
            else:
                if '^user:TaxonName' in self.data.otu_dict[key].keys():
                    value = self.data.otu_dict[key]['^user:TaxonName']
                elif '^ot:ottTaxonName' in self.data.otu_dict[key].keys():
                    value = self.data.otu_dict[key]['^ot:ottTaxonName']
                    if value is None:
                        value = self.data.otu_dict[key]['^ncbi:taxon']

            # debug("populate dict with values")
            value = str(value).replace(" ", "_")
            if value in self.sp_d:
                self.sp_d[value].append(self.data.otu_dict[key])
                # assert value = self.data.otu_dict[key]['^ncbi:taxon'] or
            else:
                self.sp_d[value] = [self.data.otu_dict[key]]
        return self.sp_d


    def make_sp_seq_dict(self, treshold, selectby):
        '''makes dict with spname as key1, key2 is gi and value is seq.
        This is used to select for the filtering step, where it
        selects how many sequences per species to keep in the alignment.
        '''
        # Note: has test, test_sp_seq_d.py, runs
        debug("in make_sp_seq_dict")
        for key in self.sp_d:
            #loop to populate dict. key1 = sp name, key2= gi number, value = seq,
            #number of items in key2 will be filtered according to threshold and already present seq
            debug(key)
            seq_d = {}
            tres_minimizer = 0
            for gi_id in self.sp_d[key]:
                if '^physcraper:status' in gi_id:
                    self.not_added = ['deleted', 'subsequence,', 'not']
                    if gi_id['^physcraper:status'].split(' ')[0] not in self.not_added:
                        if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                            if '^user:TaxonName'  in gi_id:
                                #generate entry for already existing sp
                                tres_minimizer += 1
                                user_name = gi_id['^user:TaxonName']
                                for user_name_aln, seq in self.data.aln.items():
                                    if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
                                        #debug("initial user input seq")
                                        if user_name == self.data.otu_dict[user_name_aln.label]['^user:TaxonName']:
                                            seq = seq.symbols_as_string().replace("-", "")
                                            seq = seq.replace("?", "")
                                            seq_d[user_name_aln.label] = seq
                            elif '^ot:ottTaxonName'  in gi_id:
                                #generate entry for already existing sp
                                tres_minimizer += 1
                                #debug("added in earlier round, get aln directly")
                                user_name = gi_id['^ot:ottTaxonName']
                                for user_name_aln, seq in self.data.aln.items():
                                    if '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
                                        if user_name == self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']:
                                            seq = seq.symbols_as_string().replace("-", "")
                                            seq = seq.replace("?", "")
                                            seq_d[user_name_aln.label] = seq
                        else:
                            if "^ncbi:gi" in gi_id: #this should not be needed: all new blast seq have gi
                                #debug("gi in gi_id")
                                gi_num = int(gi_id['^ncbi:gi'])
                                if gi_num in self.new_seqs.keys():
                                    seq = self.new_seqs[gi_num]
                                    seq = seq.replace("-", "")
                                    seq = seq.replace("?", "")
                                    seq_d[gi_num] = seq
                        self.sp_seq_d[key] = seq_d
        # debug(self.sp_seq_d)
        return


    def run_local_blast(self, blast_seq, blast_db):
        """run local blast to select number of sequences to be kept
        """
        # Note: has test, runs -> test_run_local_blast.py
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        out_fn = "{}_tobeblasted".format(str(blast_seq))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(blast_seq)
        os.system(cmd1)
        cmd2 = "blastn -query {} -db {}_db -out output_{}.xml -outfmt 5".format(out_fn, blast_db, out_fn)
        os.system(cmd2)
        os.chdir(general_wd)


    def calculate_mean_sd(self, hsp_scores):
        """calculates standard deviation, mean of scores which are
        used to select which sequences are selected from in the localblast
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
        mean_hsp_bits = float(bit_sum/total_seq)
        mean_sd = {"mean": mean_hsp_bits, "sd" : bit_sd}
        return mean_sd


    def read_local_blast(self, seq_d, fn):
        """reads the files of the local blast run
        and return sequenced below a threshold (0.003) and
        which are within the std of the mean scores at the moment
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
                # filter by e-value
                ## !!! maybe don"t do that....
                # if hsp.expect < 0.003:
                    gi_id = alignment.title.split(" ")[1]
                    try:
                        gi_id = int(gi_id)
                    except:
                        gi_id = gi_id
                    # debug(gi_id)
                    hsp_scores[gi_id] = {"hsp.bits" : hsp.bits, "hsp.score" : hsp.score, "alignment.length" : alignment.length, "hsp.expect" : hsp.expect}
                # else:
                #     debug("sequences highly different")
                #     debug(hsp)
        # make values to select for blast search, calculate standard deviation,mean
        mean_sd = self.calculate_mean_sd(hsp_scores)
        #select which sequences to use
        seq_blast_score = {}
        for gi_id in hsp_scores: # use only seq that are similar to mean plus minus sd
            if (hsp_scores[gi_id]["hsp.bits"] >= mean_sd["mean"]-mean_sd["sd"]) & (hsp_scores[gi_id]["hsp.bits"] <= mean_sd["mean"]+mean_sd["sd"]):
                if gi_id in seq_d:
                    seq_blast_score[gi_id] = seq_d[gi_id]
        return seq_blast_score

    def select_seq_by_local_blast(self, seq_d, fn, treshold, count):
        """select number of sequences from local blast
        to fill up to the threshold.
        only species included which have a blast score of mean plus/minus sd.
        Is used after read_local_blast
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
        # debug(self.filtered_seq)
        return self.filtered_seq

    def select_seq_by_length(self, taxon_id, treshold, count):
        """select new sequences by length"""
        debug("select_seq_by_length")

        ReqMaxLength = max(self.sp_seq_d[taxon_id].values())
        ### !!! sometimes the only seq in seq_w_maxlen is the original seq,
        # then this is the one to be added, but it will be removed,
        # later as it is no new seq! thus no new seq for that species is added
        ##

        #code is untested, tested code is the code below in the comment
        seq_w_maxlen = {}
        for key, val in self.sp_seq_d[taxon_id].iteritems():
            if len(val) == len(ReqMaxLength):
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
            subdict = {k:v for k, v in self.sp_seq_d[taxon_id].iteritems() if k not in keymax}
            second_len = max(subdict.values())
            seq2len = {}
            for key, val in subdict.iteritems():
                if len(val) == len(second_len):
                    seq2len[key] = val
            random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
            random_seq_ofsp = dict(random_seq_ofsp)
            random_seq_ofsp.update(seq_w_maxlen)
        if random_seq_ofsp != None:
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
        #         ReqMaxLength = max(self.sp_seq_d[key].values())
        #         ### !!! sometimes the only seq in seq_blast_score is the
        #         original seq, then this is the one to be added, but it will be removed,
        #         # later as it is no new seq! thus no new seq for that species is added
        #         ##
        #         seq_w_maxlen = {}
        #         for k, v in self.sp_seq_d[key].iteritems():
        #             if len(v) == len(ReqMaxLength):
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
        """add all seq to filtered_dict as the number of seq < threshold
        """
        # Note: has test, test_add_all.py: runs
        debug('add_all')
        for gi_id in self.sp_d[key]:
            # if gi_id['^physcraper:last_blasted'] != '1800/01/01':
            #     """seq was already blasted"""
            #     if '^user:TaxonName' in gi_id:
            #         user_name = gi_id['^user:TaxonName']
            #         for user_name_aln, seq in self.data.aln.items():
            #             if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
            #                 if user_name == self.data.otu_dict[user_name_aln.label]['^user:TaxonName']:
            #                     self.filtered_seq[user_name_aln.label] = seq.symbols_as_string()
            #     elif '^ot:ottTaxonName' in gi_id:
            #         user_name = gi_id['^ot:ottTaxonName']
            #         for user_name_aln, seq in self.data.aln.items():
            #             if '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
            #                 if user_name == self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']:
            #                     self.filtered_seq[user_name] = seq.symbols_as_string()
            # elif gi_id['^physcraper:last_blasted'] == '1800/01/01':
            #debug(gi_id)

            


            if '^physcraper:status' in gi_id:
                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if gi_id['^physcraper:last_blasted'] == '1800/01/01':
                        gi_num = gi_id['^ncbi:gi']
                        debug(gi_num)


                        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)


                        avg_seqlen = sum(self.data.orig_seqlen)/len(self.data.orig_seqlen) #HMMMMMMMM
                        assert self.config.seq_len_perc <= 1
                        seq_len_cutoff = avg_seqlen*self.config.seq_len_perc

                        for gi, seq in self.new_seqs.items():
                            if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                                otu_id = self.data.add_otu(gi, self.ids)
                                self.seq_dict_build(seq, otu_id, tmp_dict)


                                for gi_num, seq in tmp_dict.items():


                                    self.filtered_seq[gi_num] = seq

                        # if gi_num, seq in self.new_seqs.items():
                        #     # debug(self.new_seqs[gi_num])
                            
                        #     # seq = self.sp_seq_d[key][gi_num]
                        #     self.filtered_seq[gi_num] = seq
        return self.filtered_seq

    def loop_for_write_blast_files(self, key, selectby):
        """this loop is needed to be able to write the local blast files correctly
        """
        # Note: has test,test_loop_for_blast.py: runs
        debug("length of sp_d key")
        debug(len(self.sp_d[key]))
        nametoreturn = None
        blastfile_taxon_names = {}

        ##loop needs to happen before the other one, as we need nametoreturn in second:
        for gi_id in self.sp_d[key]:
            if '^physcraper:status' in gi_id:
                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                        if '^user:TaxonName' in gi_id:
                            user_name = gi_id['^user:TaxonName']
                            for user_name_aln, seq in self.data.aln.items():
                                if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
                                    if user_name == self.data.otu_dict[user_name_aln.label]['^user:TaxonName']:
                                        nametoreturn = user_name_aln.label
                        elif '^ot:ottTaxonName' in gi_id:
                            user_name = gi_id['^ot:ottTaxonName']
                            for user_name_aln, seq in self.data.aln.items():
                                if '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
                                    if user_name == self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']:
                                        nametoreturn = user_name_aln.label
                                        


        for gi_id in self.sp_d[key]:
            debug("in writing file for-loop")
            if '^physcraper:status' in gi_id:
                if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    # debug("giId not in not_to_add")
                    debug("genrate files used for blast")
                    if gi_id['^physcraper:last_blasted'] != '1800/01/01':
                        if '^user:TaxonName' in gi_id:
                            user_name = gi_id['^user:TaxonName']
                            for user_name_aln, seq in self.data.aln.items():
                                if '^user:TaxonName' in self.data.otu_dict[user_name_aln.label]:
                                    if user_name == self.data.otu_dict[user_name_aln.label]['^user:TaxonName']:
                                        if selectby == "blast":
                                            #debug("user")
                                            # debug(user_name)
                                            # debug(user_name_aln.label)
                                            fn = user_name_aln.label
                                            if self.downtorank != None:
                                                fn = key
                                            self.write_blast_files(fn, seq)
                                            blastfile_taxon_names[user_name] = user_name_aln.label
                                            # nametoreturn = user_name_aln.label
                        elif '^ot:ottTaxonName' in gi_id:
                            # debug("seq already in aln")
                            user_name = gi_id['^ot:ottTaxonName']
                            for user_name_aln, seq in self.data.aln.items():
                                if '^ot:ottTaxonName' in self.data.otu_dict[user_name_aln.label]:
                                    if user_name == self.data.otu_dict[user_name_aln.label]['^ot:ottTaxonName']:
                                        if selectby == "blast": # should be obsolete now?!
                                            #debug("ott")
                                            # debug(user_name_aln.label)
                                            # debug(seq)
                                            fn = user_name_aln.label
                                            if self.downtorank != None:
                                                fn = key
                                            self.write_blast_files(fn, seq)
                                            blastfile_taxon_names[user_name] = user_name_aln.label
                                            # nametoreturn = user_name_aln.label
                    else:
                        debug("make gilist as local blast database")
                        if "^ncbi:gi" in gi_id:
                            gi_num = int(gi_id['^ncbi:gi'])
                            debug(gi_num)
                            if selectby == "blast": # should be obsolete now?!
                                debug("new")
                                file_present = False
                                #debug(gi_num)
                                if gi_num in self.new_seqs.keys():
                                # for k ,v in self.new_seqs.iteritems():
                                #     if '^ncbi:gi' in v:
                                    debug("gi exists in new_seqs")
                                    # debug(v['^ncbi:gi'])
                                    # if gi_num == v['^ncbi:gi']:
                                    file_present = True
                                if file_present: # short for if file_present == True
                                    if '^physcraper:status' in gi_id:
                                        if gi_id['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                                            blast_fn = gi_num
                                            debug("write seq to db")
                                            debug(nametoreturn)
                                            seq = self.sp_seq_d[key][gi_num]
                                            if self.downtorank != None:
                                                blast_fn = key
                                                nametoreturn = key
                                            debug(blast_fn)
                                            self.write_blast_files(blast_fn, seq, db=True, fn=nametoreturn)
                                            blastfile_taxon_names[gi_num] = gi_num
                        namegi = key
        if self.downtorank != None:
            nametoreturn = key
        if nametoreturn is None:
            nametoreturn = namegi
        return nametoreturn

    def count_num_seq(self, taxon_id):
        """make some calculation of how many species are already present
        and how many seq have been found by blast for same sp. taxon concept.
        This will be used for how_many_sp_to_keep.
        """
        #this is used to exclude stuff, that has more info in sp_d,
        # but has no new seq in sp_seq_dic
        seq_present = 0
        if taxon_id in self.sp_seq_d.keys():
            for sp_keys in self.sp_seq_d[taxon_id].keys():
                if isinstance(sp_keys, str):
                    seq_present += 1
        #this calculates how many seq of species have already
        # been present in the aln
        new_taxon = True
        query_count = 0
        for item in self.sp_d[taxon_id]:
            if '^physcraper:status' in item:
                if item['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if item['^physcraper:last_blasted'] != '1800/01/01':
                        new_taxon = False
                if  item['^physcraper:status'] == "query":
                    query_count += 1
        count_dict = {'seq_present': seq_present, 'query_count': query_count, 'new_taxon': new_taxon}
        return count_dict


    def how_many_sp_to_keep(self, treshold, selectby):
        """uses the sp_seq_d and places the number of sequences
        according to threshold into the filterdseq_dict
        """
        debug("how_many_sp_to_keep")
        for taxon_id in self.sp_d:
            if len(self.sp_d[taxon_id]) <= treshold:
                ## add all stuff to self.filtered_seq[gi_n],
                # because the sequence number is below threshold,
                # thus everythng will be added
                self.add_all(taxon_id)
            elif len(self.sp_d[taxon_id]) > treshold:
                debug("filter number of sequences")
                # debug(self.sp_seq_d[taxon_id].keys())

                count_dict = self.count_num_seq(taxon_id)
                # #############
                # count = 0
                # if taxon_id in self.sp_seq_d.keys(): #if added to exclude stuff, that has more info in sp_d, but has no new seq in sp_seq_dic
                #     ## do some counting to redirect the different filters.
                #     for sp_keys in self.sp_seq_d[taxon_id].keys():
                #         # debug(type(sp_keys))
                #         if isinstance(sp_keys, str):
                #             count += 1
                #     #this calculates how many seq of species have already
                #     # been present in the aln
                #     new_taxon = True
                #     toadd = 0
                #     query_count = 0
                #     for item in self.sp_d[taxon_id]:
                #             if '^physcraper:status' in item:
                #                 not_to_add = ['deleted', 'subsequence,', 'not']
                #                 if item['^physcraper:status'].split(' ')[0] not in not_to_add:
                #                     # debug("giId not in not_to_add")
                #                     toadd += 1
                #                     if item['^physcraper:last_blasted'] != '1800/01/01':
                #                         debug("taxon not new")
                #                         new_taxon = False
                #                 if  item['^physcraper:status'] == "query":
                #                     query_count += 1
                #                     #####################


                if taxon_id in self.sp_seq_d.keys():
                    seq_present = count_dict["seq_present"]
                    query_count = count_dict["query_count"]
                    if selectby == "length":
                        self.select_seq_by_length(self.sp_seq_d[taxon_id], treshold, seq_present)
                    elif selectby == "blast":
                        # debug("blast locally")
                        # debug(seq_present)
                        # debug(treshold)
                        # debug(len(self.sp_seq_d[taxon_id]))
                        # debug("query_count")
                        # debug(query_count)
                        if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
                            debug("seq_present>0")
                            if query_count + seq_present > treshold:
                                taxonfn = self.loop_for_write_blast_files(taxon_id, selectby)
                                #species is not new in alignment, make blast with existing seq
                                # next line should be redundant
                                #if taxon_id in self.sp_d.keys():
                                for element in self.sp_d[taxon_id]:
                                    if '^ot:ottTaxonName' in element:
                                        blast_seq = "{}".format(element['^ot:ottTaxonName'])
                                        blast_seq = blast_seq.replace(" ", "_")
                                        blast_db = "{}".format(element['^ot:ottTaxonName'])
                                        blast_db = blast_db.replace(" ", "_")
                                if self.downtorank != None:
                                    taxonfn = taxon_id
                                self.run_local_blast(taxonfn, taxonfn)
                                # blast_seq = str(blast_seq.replace(" ", "_"))
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], taxonfn, treshold, seq_present)
                            elif query_count + seq_present <= treshold:
                                self.add_all(taxon_id)
                        elif seq_present == 0 and count_dict["new_taxon"] == True:
                            debug("completely new taxon to blast")
                            #species is completely new in alignment, \
                            #make blast with random species
                            # debug(self.sp_seq_d[taxon_id])
                            for item in self.sp_d[taxon_id]:
                                if '^ncbi:gi' in item:
                                    self.data.add_otu(item['^ncbi:gi'], self.ids, self.config.email)
                            blast_seq = self.sp_seq_d[taxon_id].keys()[0]
                            if type(blast_seq) == int:
                                str_db = str(taxon_id)
                            else:
                                str_db = str(blast_seq)
                            blast_db = self.sp_seq_d[taxon_id].keys()[1:]
                            # write files for local blast first:
                            seq = self.sp_seq_d[taxon_id][blast_seq]
                            self.write_blast_files(str_db, seq) #blast qguy
                            debug("blast db new")
                            debug(blast_db)
                            for blast_key in blast_db:
                                seq = self.sp_seq_d[taxon_id][blast_key]
                                self.write_blast_files(blast_key, seq, db=True, fn=str_db) #local db
                            # make local blast of sequences
                            if self.downtorank != None:
                                str_db = taxon_id
                            self.run_local_blast(str_db, str_db)
                            if len(self.sp_seq_d[taxon_id]) + seq_present >= treshold:
                                self.select_seq_by_local_blast(self.sp_seq_d[taxon_id], str_db, treshold, seq_present)
                            elif len(self.sp_seq_d[taxon_id]) + seq_present < treshold:
                                self.add_all(taxon_id)
        return


    def write_blast_files(self, user_name, seq, db=False, fn=None):
        """write local blast files which will be read by run_local_blast
        """
        debug("writing files")
        if not os.path.exists("{}/blast".format(self.data.workdir)):
            os.makedirs("{}/blast/".format(self.data.workdir))
        if db:
            fnw = "{}/blast/{}_db".format(self.workdir, fn)
            fi_o = open(fnw, 'a')
        else:
            fnw = "{}/blast/{}_tobeblasted".format(self.workdir, user_name)
            fi_o = open(fnw, 'w')
        #debug(fnw)
        fi_o.write(">{}\n".format(user_name))
        fi_o.write("{}\n".format(seq))
        fi_o.close()

    def replace_new_seq(self):
        '''replaces the self.new_seqs with the filtered_seq information
        '''
        debug("replace new seq")
        #debug(self.filtered_seq)
        keylist = self.filtered_seq.keys()
        #debug(keylist)
        keylist = [x for x in keylist if type(x) == int]
        seq_not_added = self.new_seqs.keys()
        seq_not_added = [x for x  in seq_not_added if type(x) == int]
        for gi_num in seq_not_added:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if  self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        # if self.data.otu_dict[key]['^physcraper:last_blasted'] == "1800/01/01":
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'not added, because there are enough seq per sp in tre'
        for gi_num in keylist:
            for key in self.data.otu_dict.keys():
                # debug(self.data.otu_dict[key])
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if  self.data.otu_dict[key]['^ncbi:gi'] == gi_num:
                        # if self.data.otu_dict[key]['^physcraper:last_blasted'] == "1800/01/01":
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added, as one of the representatives of the taxon'
#             self.data.otu_dict[otu_id]['^ncbi:gi'] = gi_num
        ## might not be necessary, I was missing some gi's
        # when I was replacing the original one.
        # i leave the code here for now
        # reduced_gi_dict = {k: self.data.gi_dict[k] for k in keylist}
        # debug(reduced_gi_dict)
        # self.data.gi_dict.clear()
        # self.data.gi_dict = reduced_gi_dict # data.gi_dict seems to only have newly blasted stuff
        reduced_new_seqs_dic = {k: self.filtered_seq[k] for k in keylist}
        # debug(reduced_new_seqs_dic)
        self.new_seqs = deepcopy(reduced_new_seqs_dic)
        ### !!! key is not exactly same format as before in new_seqs_otu_id
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic)
        #set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return self.new_seqs


    def add_local_seq(self, path_to_local_seq, id_to_spn_addseq_json):
        """add sequences from local database, without checking
        that it fits into the alignment,
        to your phylogeny, then do normal blast"""
        #### to DO make local blast query if new seq fit
        ##NOT WORKING, project prosponed.
        debug("add_local_seq")
        self.sp_d = {}
        # get list of sequences,
        localfiles = os.listdir(path_to_local_seq)
        for index, item in enumerate(localfiles):
            debug(index)
            debug(item)
            item = str(item)
            debug(item.startswith(".~"))
            if item.startswith(".~"):
                debug("in if")
                localfiles[index] = None
        debug(localfiles)
        localfiles = filter(None, localfiles)
        debug(localfiles)
        gi_counter = 1
        ## add assert that tests that every file is a fasta file in the folder
        for file in localfiles:
            debug(file)
            filepath = "{}/{}".format(path_to_local_seq, file)
            open_file = open(filepath)
            # for line in openFile:
            content = open_file.readlines()
            content = [x.strip() for x in content]
            content = filter(None, content) # fastest
            count = 0
            gi_list = content[::2]
            seq_list = content[1::2]
            debug(gi_list)
            debug(seq_list)
            for i in xrange(0, len(gi_list)):
                key = gi_list[i].replace(">", "")
                count = count+1
                seq = seq_list[i]
                debug(key)
                debug(seq)
                self.filtered_seq[key] = seq
                gi_counter += gi_counter
                ## numbers starting with 0000 are unpublished data
                self.data.gi_dict[key] = {'accession': "000000{}".format(gi_counter), 'title': "unpublished"} 
                debug("id_to_spn_addseq_json")
                debug(id_to_spn_addseq_json)
                self.data.otu_dict[key] = {}
                self.data.otu_dict[key]['^ncbi:gi'] = self.data.gi_dict[key]['accession']
                self.data.otu_dict[key]['^ncbi:accession'] = self.data.gi_dict[key]['accession']
                self.data.otu_dict[key]['^ncbi:title'] = self.data.gi_dict[key]['title']
                # self.otu_dict[otu_id]['^ncbi:taxon'] = ids_obj.map_gi_ncbi(gi)
                # self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ids_obj.map_gi_ncbi(gi))
                # self.otu_dict[otu_id]['^physcraper:status'] = "query"
                self.data.otu_dict[key]['^ot:ottTaxonName'] = id_to_spn_addseq_json[key]['^ot:ottTaxonName']
                self.data.otu_dict[key]['^ncbi:taxon'] = id_to_spn_addseq_json[key]['^ncbiID']
                self.data.otu_dict[key]['^ot:ottId'] = id_to_spn_addseq_json[key]['ot:ottId']
                self.data.otu_dict[key]['^physcraper:status'] = "local seq"
                self.data.otu_dict[key]['^physcraper:last_blasted'] = "1800/01/01"
                # use those seq to make a blast search, that saves one round of calculating a tree.
                debug("blast newly added seq")
                query = seq.replace("-", "").replace("?", "")
                xml_fi = "{}/{}.xml".format(self.blast_subdir, key)
                if not os.path.isfile(xml_fi):
                    sys.stdout.write("blasting seq {}\n".format(key))
                    try:
                        ###try to make a list of blast queries
                        debug(datetime.datetime.now())
                        result_handle = NCBIWWW.qblast("blastn", "nt",
                                                       query,
                                                       url_base="http://ec2-18-144-9-156.us-west-1.compute.amazonaws.com/cgi-bin/blast.cgi",
                                                       # entrez_query=equery,
                                                       hitlist_size=self.config.hitlist_size)
                        debug(datetime.datetime.now())
                        save_file = open(xml_fi, "w")
                        save_file.write(result_handle.read())
                        save_file.close()
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = str(datetime.date.today()).replace("-", "/")
                        result_handle.close()
                    except (ValueError, URLError):
                        sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(key))
        self.new_seqs = deepcopy(self.filtered_seq)
        self.new_seqs_otu_id = deepcopy(self.filtered_seq) ### !!! key is not exactly same format as before
        debug("read blast results")
        # read blast results and add to new_seqs
        for key in self.filtered_seq:
            xml_fi = "{}/{}.xml".format(self.blast_subdir, key)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    if int(alignment.title.split('|')[1]) not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                        self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                        self.data.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
                except (ValueError, URLError):
                    sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(key))
        #set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()

    def write_otu_info(self, downtorank):
        '''writes table to file with taxon names and number of representatives
        '''
        sp_d = self.sp_dict(downtorank)
        sp_info = {}
        for k in sp_d:
            sp_info[k] = len(sp_d[k])
        file = open('taxon_sampling.csv', 'w')
        file.write(str(sp_info)+"\n")
        file.close()
