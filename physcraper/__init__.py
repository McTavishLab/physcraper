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
from copy import deepcopy
from urllib2 import URLError
import configparser
import json
import pickle
from Bio.Blast import NCBIWWW, NCBIXML
import physcraper.AWSWWW as AWSWWW
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO, Entrez
from dendropy import Tree,\
                     DnaCharacterMatrix,\
                     DataSet,\
                     datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life
from peyotl.nexson_syntax import extract_tree,\
                                 extract_tree_nexson,\
                                 get_subtree_otus,\
                                 extract_otu_nexson,\
                                 PhyloSchema
from peyotl.api import APIWrapper

import collections 
import inspect
from ete2 import NCBITaxa


import datetime

def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False

_DEBUG = 1

class ConfigObj(object):
    """Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store."""
    def __init__(self, configfi):
        if _DEBUG:
             sys.stderr.write("Building config object\n")
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
        # self.blastdb = config["blast"]["localblastdb_wd"]
        self.url_base = config['blast'].get('url_base')
        if _DEBUG:
            sys.stderr.write("{}\n".format(self.email))
            sys.stderr.write("{}\n".format(self.url_base))


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
        otu_dict[otu_id]['^physcraper:last_blasted'] = "1900/01/01"
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
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir) #newick should be bare, but alignement should be DNACharacterMatrix


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
                            schema_trf="newick",
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem.
    If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.
    otu_json should encode the taxon names for each tip"""
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
   
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH

    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)

    assert(tre.taxon_namespace is aln.taxon_namespace)

    otu_newick = tre.as_string(schema=schema_trf)

    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        otu_dict = json.load(open(otu_json,"r"))
        ott_ids = [otu_dict[otu].get(u'^ot:ottId',) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)

    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)#, taxon_namespace=global_taxonnamespace)

class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match """
    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None, taxon_namespace=None ):
        #TODO add assertions that inputs are correct type!!!
        self.otu_taxonlabel_problem={}
        self.aln = alignment
        for key in self.aln.taxon_namespace: #MK_TODO this should get pull out into a function
            # print(key)
            # print(str(key))
            match = re.match("'n[0-9]{1,3}", str(key))
            # print(match)
            if match:
                # print('do something')
                newname=str(key)[2:]
                newname=newname[:-1]
                # print(newname)
                #key.label=newname
                #print(key.label)
            else: 
                newname=str(key) 
                newname=newname[:-1]
                newname=newname[1:]
                # print(newname) 
            self.otu_taxonlabel_problem[key.label]= newname
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
        self._reconcile_names()
        self.workdir = workdir #TODO - is this where the workdir should live?
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
        align with those found in the taxonomy. Spaces vs underscores
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
                self.otu_dict[tax.label]['physcraper:status'] = "deleted in name reconciliation"
            except:
                self.otu_dict[self.otu_taxonlabel_problem[tax.label]]['physcraper:status'] = "deleted in name reconciliation"
    
            self.aln.taxon_namespace.remove_taxon(tax)
        assert(self.aln.taxon_namespace==self.tre.taxon_namespace)
    

    
    def prune_short(self, min_seqlen=0):
        """Sometimes in the de-concatenating of the original alignment
        taxa with no sequence are generated.
        This gets rid of those from both the tre and the alignement. MUTATOR"""
        prune = []
        for tax, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) <= min_seqlen:
                prune.append(tax)
        if prune:
            self.aln.remove_sequences(prune)
            self.tre.prune_taxa(prune)
            #self.aln.taxon_namespace.remove_taxon(tax)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short step due to sequence shorter than {}\n".format(min_seqlen))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['physcraper:status'] = "deleted in prune short"
            self.aln.taxon_namespace.remove_taxon(tax)
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        self.reconcile()

    def reconcile(self, seq_len_perc=0.75):
        """all missing data seqs are sneaking in, but from where?!"""
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        prune = []
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        seq_len_cutoff = avg_seqlen*seq_len_perc
        for tax, seq in self.aln.items():
            # print(tax)
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(tax)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in reconcilation step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict.get(tax.label).get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['physcraper:status'] = "deleted in reconcile"
            self.remove_taxa_aln_tre(tax.label)
        aln_ids = set()
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
                self.otu_dict[leaf.taxon.label]['physcraper:status'] = "deleted due to presence in tree but not aln. ?!"
                orphaned_leafs.add(leaf)
                print(self.otu_taxonlabel_problem.keys())
                self.otu_dict[leaf.taxon.label]['physcraper:status'] = "deleted due to presence in tree but not aln. ?!"
                #orphaned_leafs.add(leaf)
               
            # else:
            #     treed_taxa.add(leaf.taxon.label)
        
        assert treed_taxa.issubset(aln_ids)
       # for key in  self.otu_dict.keys():
       #      if key not in aln_ids:
       #           sys.stderr.write("{} was in otu dict but not alignment. it should be in new seqs...\n".format(key)
        self.trim()
        self._reconciled = 1
    def trim(self, taxon_missingness = 0.75):
        '''cuts off ends of alignment, maintaining similar to original seq len
        Important bc other while whole chromosomes get dragged in!'''
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
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi]['title']
        self.otu_dict[otu_id]['^ncbi:taxon'] = ids_obj.map_gi_ncbi(gi)
        self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ids_obj.map_gi_ncbi(gi))
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = ids_obj.ott_to_name.get(self.otu_dict[otu_id]['^ot:ottId'])
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1900/01/01"#TODO check propagation...
        if _DEBUG:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi, otu_id))
        return otu_id
    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """Papara is finicky about trees and needs phylip, this writes out needed files for papara
        (except query sequences)"""
        # print("write papara files")
        # print(self.tre.taxon_namespace)
        # print(self.tre.as_string(schema="newick"))
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")#Papara is diffffffficult about root
        fi = open("{}/{}".format(self.workdir, treefilename), "w")
        fi.write(tmptre)
        fi.close()
        self.aln.write(path="{}/{}".format(self.workdir, alnfilename), schema="phylip")
    
    def write_files(self, treepath="physcraper.tre", treeschema="newick", alnpath="physcraper.fas", alnschema="fasta"):
        """Outputs both the streaming files and a ditechecked"""
        #First write rich annotation json file with everything needed for later?
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema,
                       unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)
    def write_otus(self, filename, schema='table'):
        """Writes out OTU dict as json"""
        assert schema in ['table', 'json']
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)
    def remove_taxa_aln_tre(self, taxon_label):
        #self.aln.taxon_namespace=self.tre.taxon_namespace
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        #assert(self.aln.taxon_namespace==self.tre.taxon_namespace)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.taxon_namespace.remove_taxon(tax)
            self.tre.prune_taxa([tax])
            self.otu_dict[tax.label]['physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['physcraper:status'] = "deleted, but it wasn't in teh alignemnet..."
        # print(self.tre)
        # print(self.tre.taxon_namespace)
    

    def write_labelled(self, label, treepath="labelled.tre", alnpath="labelled.fas"):
        """output tree and alignement with human readable labels
        Jumps through abunch of hoops to make labels unique.
        NOT MEMORY EFFICIENT AT ALL"""
        # print("write labelled files")
        assert label in ['^ot:ottTaxonName','user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
        # print(label)
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        tmp_fasta = self.aln.as_string(schema="fasta")

        tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta,
                                         schema="fasta")
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:      
            new_label = self.otu_dict[taxon.label].get(label)
            if new_label:
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
                new_names.add(new_label)
                taxon.label = new_label
            elif self.otu_dict[taxon.label].get("user:TaxonName"):
                new_label = self.otu_dict[taxon.label].get("user:TaxonName")
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
                new_names.add(new_label)
                taxon.label = new_label    
            elif self.otu_dict[taxon.label].get("^ot:originalLabel"):
                new_label = self.otu_dict[taxon.label].get("^ot:originalLabel")
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
                new_names.add(new_label)
                taxon.label = new_label
            elif self.otu_dict[taxon.label].get("^ncbi:taxon"):
                new_label = " ".join(["ncbi", str(self.otu_dict[taxon.label].get("^ncbi:taxon"))])
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
    def dump(self, filename = "att_checkpoint.p"):
#        frozen = jsonpickle.encode(self)
#        with open('{}/{}'.format(self.workdir, filename), 'w') as pjson:
#            pjson.write(frozen)
        pickle.dump(self, open("{}/{}".format(self.workdir,filename), "wb" ) )
        #TODO... write as proper nexml?!

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


class IdDicts(object):
    """Wraps up the annoying conversions"""#TODO - could - should be shared acrosss runs?! .... nooo.
    def __init__(self, config_obj, workdir):
        """Generates a series of name disambiguation dicts"""
        self.workdir = workdir
        self.config = config_obj
        self.ott_to_ncbi = {}
        self.ncbi_to_ott = {}
        self.ott_to_name = {}
        self.gi_ncbi_dict = {}
        fi = open(config_obj.ott_ncbi) #TODO need to keep updated
        for lin in fi:
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
    def map_gi_ncbi(self, gi):
        """get the ncbi taxon id's for a gi input"""
# <<<<<<< HEAD
# #        mapped_taxon_ids = open("{}/id_map.txt".format(self.workdir), "a")
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
#             if _DEBUG:
#                 print("get taxid")
#             Entrez.email = self.config.email
#             # print(Entrez.email)
#             # print(gi)
#             handle = Entrez.efetch(db="nucleotide",id=gi, retmode="xml")
#             tax_name = Entrez.read(handle)[0]['GBSeq_source']

# ### !!!! NCBITaxa seems to be a bit slower than the genbank query but more reliabale!
#             ncbi = NCBITaxa()    
#             tax_info = ncbi.get_name_translator([tax_name])
#             tax_id = tax_info.items()[0][1][0]
            try:
                Entrez.email = self.config.email
                handle = Entrez.efetch(db="nucleotide",id=gi, retmode="xml")
                tax_name = Entrez.read(handle)[0]['GBSeq_source']
                ncbi = NCBITaxa()    
                tax_info = ncbi.get_name_translator([tax_name])
                tax_id = tax_info.items()[0][1][0]
                self.gi_ncbi_dict[gi] = tax_id
            except:
                if _DEBUG:
                    sys.stderr.write("Entrez fetch failed, for GI:{} running subprocess to get taxon id\n".format(gi))
                tax_id = int(subprocess.check_output(["bash", self.config.get_ncbi_taxonomy,
                                                      "{}".format(gi),
                                                      "{}".format(self.config.ncbi_dmp)]).split('\t')[1])  
                self.gi_ncbi_dict[gi] = tax_id    
                       
# =======
        # mapped_taxon_ids = open("{}/id_map.txt".format(self.workdir), "a")
        # if _DEBUG:
        #     sys.stderr.write("Mapping gi {} to taxon\n".format(gi))
        # if gi in self.gi_ncbi_dict:
        #     tax_id = int(self.gi_ncbi_dict[gi])
        # else:
        #     tax_id = None


# #############################3

# >>>>>>> aws_blast
#             # try:
#             #     Entrez.email = self.config.email
#             #     handle = Entrez.efetch(db="nucleotide",id=gi, retmode="xml")
#             #     tax_name = Entrez.read(handle)[0]['GBSeq_source']
#             #     ncbi = NCBITaxa()    
#             #     tax_info = ncbi.get_name_translator([tax_name])
#             #     tax_id = tax_info.items()[0][1][0]
#             #     self.gi_ncbi_dict[gi] = tax_id
#             # except:
#             #     if _DEBUG:
#             #         sys.stderr.write("Entrez fetch failed, for GI:{} running subprocess to get taxon id\n".format(gi))
#             #     tax_id = int(subprocess.check_output(["bash", self.config.get_ncbi_taxonomy,
#             #                                           "{}".format(gi),
#             #                                           "{}".format(self.config.ncbi_dmp)]).split('\t')[1])  
#             #     self.gi_ncbi_dict[gi] = tax_id    
#         if tax_id == None:
#             sys.stderr.write("Taxon not identified GI:{}\n".format(gi))             
#             tax_id = 0
#         if _DEBUG:
#             sys.stderr.write("gi {} mapped to taxon ncbi {}\n".format(gi, tax_id))
#         return tax_id
    def dump(self):
        filename = "id_pickle.p"
        pickle.dump(IdDicts, open("{}/{}".format(self.workdir,filename), "wb" ) )
        

class PhyscraperScrape(object): #TODO do I wantto be able to instantiate this in a different way?!
    #set up needed variables as nones here?!
    #TODO better enforce ordering
    """This is the class that does the perpetual updating"""
    def __init__(self, data_obj, ids_obj, config_obj):
        #todo check input types assert()
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = config_obj
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
                    equery = "{}:{}[mdat]".format(self.mrca_ncbi,
                                                  last_blast,
                                                  today)
                    query = seq.symbols_as_string().replace("-", "").replace("?", "")
                    
                    fn_oldseq = taxon.label +  "_tobeblasted"
                    print("{}/{}".format(self.blast_subdir,fn_oldseq))
                    fi_old = open("{}/{}".format(self.blast_subdir,fn_oldseq), 'w')
                    fi_old.write(">{}\n".format(taxon.label))
                    fi_old.write("{}\n".format(query))
                    fi_old.close()

                    xml_fi = "{}/{}.xml".format(self.blast_subdir, fn_oldseq)
                    if not os.path.isfile(xml_fi):
                        sys.stdout.write("blasting seq {}\n".format(taxon.label))
# <<<<<<< HEAD
                        try:
                            # result_handle = NCBIWWW.qblast("blastn", "nt",
                            #                                query,
                            #                                entrez_query=equery,
                            #                                hitlist_size=self.config.hitlist_size)
                            print(datetime.datetime.now())
                            blastcmd = "blastn -query " + "{}/{}".format(self.blast_subdir,fn_oldseq) +" -db /shared/localblastdb_meta/nt -out " + xml_fi + " -outfmt 5 -num_threads 8"
                            
                            # blastcmd = "blastn -query " + "{}/{}".format(self.blast_subdir,fn_oldseq) +" -db " + self.blastdb +"/nt -out " + xml_fi + " -outfmt 5 -num_threads 8"
                            
                            print(blastcmd)
                            os.system(blastcmd)
                            


                            # save_file = open(xml_fi, "w")
                            # save_file.write(result_handle.read())
                            # save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                            # result_handle.close()
#                         except (ValueError, URLError):
#                             sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
# =======
                       # try:
                        except:
                            if self.config.url_base:
                                result_handle = AWSWWW.qblast("blastn",
                                                               "nt",
                                                               query,
                                                               url_base = self.config.url_base,
                                                               entrez_query=equery,
                                                               hitlist_size=self.config.hitlist_size,
                                                               num_threads=4)
                            else:
                                result_handle = AWSWWW.qblast("blastn",
                                                               "nt",
                                                               query,
                                                               entrez_query=equery,
                                                               hitlist_size=self.config.hitlist_size)
                            save_file = open(xml_fi, "w")
                            save_file.write(result_handle.read())
                            save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                            result_handle.close()
                       # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                       #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
        self._blasted = 1
        return
    def read_blast(self, blast_dir=None):
        """reads in and prcesses the blast xml files"""
        if not blast_dir:
            blast_dir = self.blast_subdir
        if not self._blasted:
            self.run_blast()
        for taxon in self.data.aln:
            # xml_fi =  str(taxon.label) + ".xml"
            # xml_fi = "{}/{}.xml".format(blast_dir, taxon.label)
            xml_fi = "{}/{}_tobeblasted.xml".format(self.blast_subdir, taxon.label)
                    
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    # print(self.data.gi_dict)
                                    # print(int(alignment.title.split('|')[1]))
                                    # print(alignment.title.split("|"))
                                    gi = alignment.title.split("|")[1]
                                    # print(gi)
                                    # print(something_stupid)
                                    # print(blast_record.query_id)
                                    # print(alignment.hit_id)
                                    # print("hsp")
                                    # print(hsp)
                                    # print(something_stupid)
                                    if gi not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                        self.new_seqs[gi] = hsp.sbjct
                                        self.data.gi_dict[gi] = alignment.__dict__
                                    # if int(alignment.title.split('|')[1]) not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                        # self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                        # self.data.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
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
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        for tax_lab in tax_list:
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(inc_seq) >= len(new_seq):
                if inc_seq.find(new_seq) != -1:
                    sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                    self.data.otu_dict[tax_lab]['physcraper:status'] = "subsequence, not added"
                    return
            else:
                if new_seq.find(inc_seq) != -1:#
                    if self.data.otu_dict[tax_lab].get('^physcraper:status') == "original":
                        # print("delete because its a superseq")

                        sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax_lab))
                        self.data.otu_dict[tax_lab]['physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        return
                    else:
                        # print("delete because its a superseq")
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[tax_lab]['physcraper:status'] = "new seq added in place of {}".format(tax_lab)
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
        print("remove identical seqs")
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        #Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.data.orig_seqlen)/len(self.data.orig_seqlen) #HMMMMMMMM
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


    def dump(self, filename = "scrape_checkpoint.p"):
#        frozen = jsonpickle.encode(self)
#        with open('{}/{}'.format(self.workdir, filename), 'w') as pjson:
#            pjson.write(frozen)
        pickle.dump(self, open("{}/{}".format(self.workdir,filename), "wb" ) )
        #TODO... write as proper nexml?!    

    def how_many_sp_to_keep(self, treshold):
        if _DEBUG:
            sys.stderr.write("Runing func: how_many_sp_to_keep\n")
        sp_in_aln = self.sp_d
        #sp_blast = self.new_seqs
        if _DEBUG:
            sys.stderr.write("ncbi taxon id with blast seq info in list {}\n".format(sp_in_aln))
            sys.stderr.write("gi and new_seqs{}\n".format(self.newseqs))
        # print('give otu id and seq')
        # print(self.new_seqs_otu_id)
        # print(" next is dict with otu id key and all other numbers")
        # print(self.data.otu_dict)
            sys.stderr.write("self.seq_dict_build")
            sys.stderr.write(self.seq_dict_build)
            sys.stderr.write(dir(self.seq_dict_build))
        # read blast data to check similarity of sequences: copy of read_blast
            sys.stderr.write(self.seq_dict_build.__getattribute__)
        blast_dir = self.blast_subdir        
        self.new_seqs_dic = {}
        for key in self.sp_d:
            gi_l = []
            print(key)
            print(self.sp_d[key])

            print("get infos of sp_d")
            print(len(self.sp_d[key]))

            if len(self.sp_d[key]) > treshold:
                print(treshold)
                print("There are more sequences for the taxon than your treshold, will remove the shortest.")
                seq_d = {}
                print(self.sp_d[key])
                for item in self.sp_d[key]:
                    print(item)
                    if '^ncbi:gi' in item: 
                        print("make gi list")
                        gi = item['^ncbi:gi']
                        print(gi)
                        gi_l.append(gi)

                        for gi in gi_l:
                            seq = self.new_seqs[gi]
                            print(type(seq))
                            seq = seq.replace("-", "")
                            seq_d[gi] = seq

                    elif '^user:TaxonName' in item: 
                        print("no gi, get aln directly")
                        userName = item['^user:TaxonName']
                        print(userName)
                        for userName, seq in self.data.aln.items():
                            print(type(seq))
                            seq = seq.symbols_as_string().replace("-", "")
                            #seq_d[userName.label] = [seq, len(seq)]  
                            seq_d[userName.label] = seq 

                print(gi_l)
                print("in seq_d")
                print(seq_d)

                print("seq_d.values()")

                print(seq_d.values())
                maxLength = max(seq_d.values())
                print("maxLength")

                print(maxLength)
                key_for_max = [key for key in seq_d if len(seq_d[key])==len(maxLength)]
                print([v for v in seq_d.values() if len(v) == len(maxLength)])
                #seq_w_maxlen = [v for v in seq_d.values() if len(v) == len(maxLength)]
                seq_w_maxlen ={}
                for k,v in seq_d.iteritems():
                    if len(v) == len(maxLength):
                        seq_w_maxlen[k] = v 


                print(seq_w_maxlen)
                print('start loop')
                print(len(seq_w_maxlen))
                if len(seq_w_maxlen) > treshold:
                    import random
                    random_seq_ofsp = random.sample(seq_w_maxlen.items(), treshold)
                    random_seq_ofsp = dict(random_seq_ofsp)

                    # d = seq_w_maxlen.fromkeys(range(len(seq_w_maxlen)))
                    # keys = random.sample(list(d), treshold)
                    # values = [d[k] for k in keys]


                elif len(seq_w_maxlen) == treshold:
                    random_seq_ofsp = seq_w_maxlen

                else:
                    toselect = range(len(seq_w_maxlen), treshold)   
                    print("secondlen")
                    sort_seq_d_v = sorted(seq_d.values())
                    print(sort_seq_d_v)
                    secondlen = sort_seq_d_v[1]
                    print(secondlen)
                    SeqToUse = {}
                    print(range(1,treshold))
                    for val in range(1,treshold):
                        print(val)
                        random_seq_ofsp = sort_seq_d, key=seq_d.get([val])

                for key in random_seq_ofsp.keys():
                    print(key)
                    self.new_seqs_dic[key] = random_seq_ofsp[key]
            else: 
                for item in self.sp_d[key]:
                    if '^ncbi:gi' in item: 
                        gi = item['^ncbi:gi']
                        gi_l.append(gi)

                        for gi in gi_l:
                            seq = self.new_seqs[gi]
                            seq = seq.replace("-", "")
                            self.new_seqs_dic[gi] = seq

                    elif '^user:TaxonName' in item: 
                        userName = item['^user:TaxonName']
                        for userName, seq in self.data.aln.items():
                            seq = seq.symbols_as_string().replace("-", "")
                            #seq_d[userName.label] = [seq, len(seq)]  
                            self.new_seqs_dic[userName.label] = seq 

    def replace_new_seq(self):
        # print("self.data.add_otu")
        # print(self.data.add_otu)
        
        # print(something)


        print("self.data.gi_dict")
        print(self.data.gi_dict)






        keylist = self.new_seqs_dic.keys()
        keylist = [x for x in keylist if type(x) == int]
        print(keylist)

        reduced_gi_dict = {k: self.data.gi_dict[k] for k in keylist}
        print(reduced_gi_dict)
        self.data.gi_dict.clear()
        self.data.gi_dict = reduced_gi_dict
        print("before clear:")
        print(self.new_seqs)
        # print(something)
        print("random_seq_ofsp")
        self.new_seqs.clear()
        print(self.new_seqs_dic)
        print(self.new_seqs)
        for key in self.new_seqs_dic.keys():
            print(key)
            print(type(key))
            if type(key)==int: 
                self.new_seqs[key] = self.new_seqs_dic[key]



    def sp_dict(self):
        print("make sp_dict")
        self.sp_d = {}
        for key in self.data.otu_dict:
            print("key")
            print(key)
            
            if '^user:TaxonName' in self.data.otu_dict[key].keys():
                print("in if")

                value = self.data.otu_dict[key]['^user:TaxonName']
                print(value)
            elif '^ot:ottTaxonName' in self.data.otu_dict[key].keys():
                print("in elif")
                value = self.data.otu_dict[key]['^ot:ottTaxonName']
                print(value)
            elif' ^ncbi:taxon' in self.data.otu_dict[key].keys():
                print("in elif2")

                value = self.data.otu_dict[key]['^ncbi:taxon']
                print(value)
            else:   
                print("it breaks")
                print(self.data.otu_dict[key])
                

            if value in self.sp_d:    
                self.sp_d[value].append(self.data.otu_dict[key])
            else:
                self.sp_d[value] = [self.data.otu_dict[key]]
    
        print(self.sp_d) 
        #self.sp_d
        return

    def write_query_seqs(self):
        """writes out the query sequence file"""
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
        sys.stdout.write("align query seqs")

        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        sys.stdout.write("aligning query sequences \n")
        self.data.write_papara_files()
        os.chdir(self.workdir)#Clean up dir moving
        print("trying to call papara")
        try:
            print("I call papara")

            pp = subprocess.call(["papara",
                                  "-t", "random_resolve.tre",
                                  "-s", "aln_ott.phy",
                                  "-q", self.newseqs_file,
                                  "-n", papara_runname]) #FIx directory ugliness
            sys.stdout.write("Papara done")
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it insatlled?\n")
                sys.exit()
            # handle file not found error.
            else:
            # Something else went wrong while trying to run `wget`
                raise
        os.chdir('..')
        assert os.path.exists(path="{}/papara_alignment.{}".format(self.workdir, papara_runname))
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.{}".format(self.workdir, papara_runname), schema="phylip")
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
        #print(os.getcwd())
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
        os.chdir('..')
        self._query_seqs_placed = 1
    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree"""
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        p2 = subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                              "-s", "papara_alignment.extended",
                              "-t", "place_resolve.tre",
                              "-p", "1",
                              "-n", "{}".format(self.date)])
        os.chdir('..')#TODO mordir not always one down!
        self._full_tree_est = 1
    def generate_streamed_alignment(self, treshold=None):
        """runs the key steps and then replaces the tree and alignme nt with the expanded ones"""
        self.read_blast()
        self.data.dump()
	if treshold!=None:
	        self.sp_dict()
        	self.how_many_sp_to_keep(treshold)
        	self.replace_new_seq()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
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
        self.reset_markers()
        self.data.dump()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))
