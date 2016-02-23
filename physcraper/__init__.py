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
import pickle
import unicodedata
from copy import deepcopy
import configparser
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO, Entrez
from dendropy import Tree, DnaCharacterMatrix
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees, iter_otus
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree, extract_tree_nexson, get_subtree_otus, extract_otu_nexson, PhyloSchema
from peyotl.api import APIWrapper

class StudyInfo(object):
    """Wraps up tree, alignment, and study info"""
    def __init__(self,
                 study_id,
                 tree_id,
                 seqaln,
                 mattype):
        self.study_id = study_id
        self.tree_id = tree_id
        self.seqaln = seqaln
        self.mattype = mattype




class PhyscraperSetup(object):
    """This needs to vet the inputs, standardize names, and prep for
    first run through of blasting. Does not blast itself!"""
    _found_mrca = 0
    _id_dicts = 0
    _study_get = 0
    _phylesystem = 0
    def __init__(self,
                 study_info,
                 runname,
                 configfi='/home/ejmctavish/projects/otapi/physcraper/config'):
        """initialized object, most attributes generated through self._checkArgs using config file."""
        self.study_id = study_info.study_id
        self.tree_id = study_info.tree_id
        self.seqaln = study_info.seqaln
        self.mattype = study_info.mattype
        self.runname = runname
        self.configfi = configfi
        assert os.path.isfile(configfi)
        self._read_config()
        self._get_study()
        self.ott_to_ncbi = {}
        self.ott_to_name = {}
        self.ncbi_to_ott = {}
        self.phy = None
        self.gi_ncbi_dict = {}
        self.mrca_ott = None
        self.aln = None
        self.orig_aln = None
        self.tre = None
        self.mrca_ncbi = None
        self.orig_seqlen = None
        try:
            self._reconcile_names()
        except:
            print "name reconcillation error"
        self.workdir = runname
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
    def _read_config(self):
        self.config = configparser.ConfigParser()
        self.config.read(self.configfi)
        self.e_value_thresh = self.config['blast']['e_value_thresh']
        self.hitlist_size = int(self.config['blast']['hitlist_size'])
        self.seq_len_perc = float(self.config['physcraper']['seq_len_perc'])
        self.get_ncbi_taxonomy = self.config['ncbi.taxonomy']['get_ncbi_taxonomy']
        self.ncbi_dmp = self.config['ncbi.taxonomy']['ncbi_dmp']
    def _make_id_dicts(self):
        """Generates a series of name disambiguation dicts"""
        fi = open(self.config['ncbi.taxonomy']['ott_ncbi']) #TODO need to keep updated
        for lin in fi:
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()
        fi.close()
        if os.path.isfile("{}/id_map.txt".format(self.workdir)): #todo config?!
            fi = open("{}/id_map.txt".format(self.workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]
        self._id_dicts = 1
    def _get_study(self):
        """gets the nexson from phylesystem"""
        if not self._phylesystem:
            self._phylesystem_setup()
        self.nexson = self.phy.get_study(self.study_id)['data']
        #print(self.nexson)
        self._study_get = 1
    def _get_mrca(self):
        """finds the mrca of the taxa in the ingroup of the original
        tree. The blast serach later is limited to descendents of this
        mrca according to the ncbi taxonomy"""
        if not self._study_get:
            self._get_study()
        self._make_id_dicts()
        orig_ott_ids = get_subtree_otus(self.nexson,
                                        tree_id=self.tree_id,
                                        subtree_id="ingroup",
                                        return_format="ottid")
        if None in orig_ott_ids:
            orig_ott_ids.remove(None)
        mrca_node = tree_of_life.mrca(ott_ids=list(orig_ott_ids), wrap_response=True)
        self.mrca_ott = mrca_node.nearest_taxon.ott_id
        self.mrca_ncbi = self.ott_to_ncbi[self.mrca_ott]
        self._found_mrca = 1
    def _phylesystem_setup(self):
        """Setups up the phylsesystem object"""
        phylesystem_loc = self.config['phylesystem']['location']
        self.phy = PhylesystemAPI(get_from=phylesystem_loc)
        self._phylesystem = 1
    def _reconcile_names(self):
        """This checks that the tree "original labels" from phylsystem
        align with those found in the taxonomy. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!"""
        otus = get_subtree_otus(self.nexson, tree_id=self.tree_id)
        self.treed_taxa = {}
        self.otu_dict = {}
        self.orig_lab_to_otu = {}
        for otu_id in otus:
            self.otu_dict[otu_id] = extract_otu_nexson(self.nexson, otu_id)[otu_id]
            self.otu_dict[otu_id]['^physcraper:status'] = "original"
            self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1900/01/01"
            orig = self.otu_dict[otu_id].get(u'^ot:originalLabel').replace(" ", "_")
            self.orig_lab_to_otu[orig] = otu_id
            self.treed_taxa[orig] = self.otu_dict[otu_id].get(u'^ot:ottId')
        self.aln = DnaCharacterMatrix.get(path=self.seqaln, schema=self.mattype)
        for tax in self.aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")#Forcing all spaces to underscore
        missing = [i.label for i in self.aln.taxon_namespace if i.label not in self.treed_taxa.keys()]
        if missing:
            errmf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            raise ValueError(errm)
    def _prune(self):
        """Sometimes in the de-concatenating of the original alignment
        taxa with no sequence are generated.
        This gets rid of those from both the tre and the alignement"""
        prune = []
        tmp_dict = {}
        self.orig_seqlen = []
        for taxon, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) == 0:
                prune.append(taxon.label)
            else:
                tmp_dict[taxon.label] = seq
                self.orig_seqlen.append(len(seq.symbols_as_string().translate(None, "-?")))
        newick = extract_tree(self.nexson,
                              self.tree_id,
                              PhyloSchema('newick',
                                          output_nexml2json='1.2.1',
                                          content="tree",
                                          tip_label="ot:originalLabel"))
        newick = newick.replace(" ", "_") #UGH Very heavy handed
        self.orig_aln = self.aln
        self.aln = DnaCharacterMatrix.from_dict(tmp_dict)
        self.tre = Tree.get(data=newick,
                            schema="newick",
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        self.tre.prune_taxa_with_labels(prune)
        for tax in prune:
            del self.otu_dict[self.orig_lab_to_otu[tax]]
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'w')
            fi.write("taxa pruned from tree and alignment due to excessive missing data\n")
            for tax in prune:
                fi.write("{}\n".format(tax))
            fi.close()
        for taxon in self.aln.taxon_namespace:
            otu_id = self.orig_lab_to_otu[taxon.label]
            taxon.label = otu_id.encode('ascii')

    def write_files(self):
        self.tre.write(path="{}/original_pruned.tre".format(self.workdir),
                       schema="newick",
                       unquoted_underscores=True)
        self.aln.write(path="{}/original_pruned.fas".format(self.workdir),
                       schema="fasta")
    def __getstate__(self):
        """hacky way to fix pickling"""
        result = self.__dict__.copy()
        del result['config']
        del result['phy']
        return result
    def setup_physcraper(self):
        """run the key steps"""
        self._read_config()
        self._get_mrca()
        self._reconcile_names()
        self._prune()
        self.write_files()
# I think it should write out to files so that the PHyscarpe scrape object can just read the aln etc from those, maybe the shaould be where the OTU dict comes from too...
#Hmmm.

class PhyscraperScrape(object): #TODO do I wantto be able to instantiate this in a different way?!
    #set up needed variables as nones here?!
    #TODO better enforce ordering
    """This is the class that does the perpetual updating"""
    
    def __init__(self, setup_obj=None):
        if setup_obj:#this should need to happen ony the first time the scrape object is run.
            self.workdir = deepcopy(setup_obj.workdir)
            self.runname = deepcopy(setup_obj.runname)
#            if os.path.exists("{}/aln_ott.phy".format(self.workdir)):
#                self.aln = DnaCharacterMatrix.get(path="{}/aln_ott.phy".format(self.workdir), schema="phylip")
#                self.tre = Tree.get(path="{}/random_resolve.tre".format(self.workdir),
#                                    schema="newick",
#                                    preserve_underscores=True,
#                                    taxon_namespace=self.aln.taxon_namespace)
#                print "PULLING FROM FILE"
#            else:
            self.aln = deepcopy(setup_obj.aln)
            self.tre = deepcopy(setup_obj.tre) #
            self.otu_dict = deepcopy(setup_obj.otu_dict)
            self.ott_to_ncbi = deepcopy(setup_obj.ott_to_ncbi)
            self.ncbi_to_ott = deepcopy(setup_obj.ncbi_to_ott)
            self.gi_ncbi_dict = deepcopy(setup_obj.gi_ncbi_dict)
            self.ott_to_name = deepcopy(setup_obj.ott_to_name)
            self.mrca_ncbi = deepcopy(setup_obj.mrca_ncbi)
            self.e_value_thresh = deepcopy(setup_obj.e_value_thresh)
            self.hitlist_size = deepcopy(setup_obj.hitlist_size)
            self.seq_len_perc = deepcopy(setup_obj.seq_len_perc)
            self.get_ncbi_taxonomy = deepcopy(setup_obj.get_ncbi_taxonomy)
            self.ncbi_dmp = deepcopy(setup_obj.ncbi_dmp)
            self.orig_seqlen = deepcopy(setup_obj.orig_seqlen)
            self.ps_otu = 1
            self.log = "{}/logfile".format(self.workdir)
        self.new_seqs = {}
        self.new_seqs_otu_id = {}
        self.otu_by_gi = {}
        self.gi_dict = {}
        self._to_be_pruned = []
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        if os.path.exists(self.tmpfi):
            dat = open(self.tmpfi).readline() #finsh what was started!
            self.today = dat.strip()
        else:
            self.today = str(datetime.date.today())
            fi = open(self.tmpfi, 'w')
            fi.write(self.today)
            fi.close()
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if os.path.exists("{}/last_completed_update".format(self.workdir)):
            last = open("{}/last_completed_update".format(self.workdir)).readline()
            self.lastupdate = last.strip()
        else:
            self.lastupdate = '1900/01/01'
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "{}/{}.fasta".format(self.workdir, self.today)
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
    def _write_files(self):
        """Outputs both the streaming files and a ditechecked"""
        #First write rich annotation json file with everything needed for later?
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;",";")#Papara is diffffffficult about root
        fi = open("{}/random_resolve.tre".format(self.workdir), "w")
        fi.write(tmptre)
        fi.close()
        self.aln.write(path="{}/aln_ott.phy".format(self.workdir),
                       schema="phylip")
        self.aln.write(path="{}/aln_ott.fas".format(self.workdir),
                       schema="fasta")
        self.tre.write(path="{}/random_resolve{}.tre".format(self.workdir, self.today),
                       schema="newick",
                       unquoted_underscores=True)
        self.aln.write(path="{}/aln_ott{}.fas".format(self.workdir, self.today),
                       schema="fasta")
    def write_labelled(self, label='^ot:ottTaxonName'):
        """output tree and alignement with human readble labels"""
        assert label in ['^ot:ottTaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]

        tmp_tre = Tree.get(path="{}/random_resolve{}.tre".format(self.workdir, self.today),
                            schema="newick",
                            preserve_underscores=True)
        tmp_aln = DnaCharacterMatrix.get(path="{}/aln_ott.fas".format(self.workdir), 
                                         schema="fasta",
                                         taxon_namespace=tmp_tre.taxon_namespace)
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            new_label = self.otu_dict[taxon.label].get(label)
            if new_label:
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
        tmp_tre.write(path="{}/labelled{}.tre".format(self.workdir, self.today),
                      schema="newick", 
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path="{}/labelled{}.fas".format(self.workdir, self.today),
                      schema="fasta")
    def run_blast(self): #TODO Should this be happening elsewhere?
        '''generates the blast queries and sames them to xml'''
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
#        equery = "txid{}[orgn]".format(self.mrca_ncbi)
        for taxon, seq in self.aln.items():
            otu_id = taxon.label
            last_blast = self.otu_dict[otu_id]['^physcraper:last_blasted']
            today = str(datetime.date.today()).replace("-", "/")
            #sys.stdout.write("last blast {}\n".format(last_blast))
            #sys.stdout.write("today {}\n".format(today))
            if abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days) > 14: #TODO make configurable
                #sys.stdout.write("date is {}\n".format(abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days)))
                equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi,
                                                               last_blast,
                                                               today)
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
                if not os.path.isfile(xml_fi):
                    result_handle = NCBIWWW.qblast("blastn", "nt",
                                                   query,
                                                   entrez_query=equery,
                                                   hitlist_size=self.hitlist_size)
                    save_file = open(xml_fi, "w")
                    save_file.write(result_handle.read())
                    save_file.close()
                    self.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                    result_handle.close()
        self._blasted = 1
#        self.lastupdate = self.today
        return #TODO pull the repeated runname
    def read_blast(self, blast_dir=None):
        '''reads in and prcesses the blast xml files'''
        if not blast_dir:
            blast_dir = self.blast_subdir
        if not self._blasted:
            self.run_blast()
        for taxon in self.aln:
            xml_fi = "{}/{}.xml".format(blast_dir, taxon.label)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                blast_records = NCBIXML.parse(result_handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if float(hsp.expect) < float(self.e_value_thresh):
                                if int(alignment.title.split('|')[1]) not in self.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                    self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                    self.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
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
        for tax in tax_list:
            inc_seq = seq_dict[tax].replace("-", "")
            if len(inc_seq) >= len(new_seq):
                if inc_seq.find(new_seq) != -1:
                    sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax))
                    del self.otu_dict[label]
                    return
            else:
                if new_seq.find(inc_seq) != -1:#
                    if self.otu_dict[tax].get('^physcraper:status') == "original":
                        sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax))
                        seq_dict[label] = seq
                        return
                    else:
                        del seq_dict[tax]
                        self.tre.prune_taxa_with_labels(tax)
                        seq_dict[label] = seq
                        del self.otu_dict[tax]  #hmmmmmmmm?! CHECK ME!!
#                        print "TAXON {} should no longer be in OTU dict!"
                        sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax, label, tax))
                        return
        sys.stdout.write(".")
        seq_dict[label] = seq
        return
    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        Does not test if they are identical to ones in the original alignment...."""
        tmp_dict = dict((taxon.label, self.aln[taxon].symbols_as_string()) for taxon in self.aln)
        old_seqs = tmp_dict.keys()
        #Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        seq_len_cutoff = avg_seqlen*self.seq_len_perc
        for gi, seq in self.new_seqs.items():
            if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                otu_id = self._add_otu(gi)
                self.seq_dict_build(seq, otu_id, tmp_dict)
        for tax in old_seqs:
            try:
                del tmp_dict[tax]
            except KeyError:
                pass
        self.new_seqs_otu_id = tmp_dict # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        sys.stdout.write("{} new sequences added from genbank\n".format(len(self.new_seqs)))
    def _add_otu(self, gi):
        """generates an otu_id for new sequences and adds them into the otu_dict"""
        otu_id = "otuPS{}".format(self.ps_otu)
        self.ps_otu += 1
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi]['title']
        self.otu_dict[otu_id]['^ncbi:taxon'] = self.map_gi_ncbi(gi)
        self.otu_dict[otu_id]['^ot:ottId'] = self.ncbi_to_ott.get(self.map_gi_ncbi(gi))
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = self.ott_to_name.get(self.otu_dict[otu_id]['^ot:ottId'])
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1900/01/01"
        return otu_id
    def map_gi_ncbi(self, gi):
        """get the ncbi taxon id's for a gi input"""
        mapped_taxon_ids = open("{}/id_map.txt".format(self.workdir), "a")
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
            try:
                tax_id = int(subprocess.check_output(["bash", self.get_ncbi_taxonomy,
                                                      "{}".format(gi),
                                                      "{}".format(self.ncbi_dmp)]).split('\t')[1])
            except ValueError:
                sys.stderr.write("ncbi_dmp file needs to be updated. Do so and rerun.")
            mapped_taxon_ids.write("{}, {}\n".format(gi, tax_id))
            self.gi_ncbi_dict[gi] = tax_id
            assert tax_id  #if this doesn't work then the gi to taxon mapping needs to be updated - shouldhappen anyhow perhaps?!
        mapped_taxon_ids.close()
        return tax_id
    def write_query_seqs(self):
        """writes out the query sequence file"""
        if not self._blast_read:
            self.read_blast()
        fi = open(self.newseqs_file, 'w')
        sys.stdout.write("writing out sequences\n")
        for otu_id in self.new_seqs_otu_id.keys():
            if otu_id not in self.aln: #new seqs only
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1
    def align_query_seqs(self, papara_runname="extended"):
        """runs papara on the tree, the alinment and the new query sequences"""
        if not self._query_seqs_written:
            self.write_query_seqs()
        sys.stdout.write("aligning query sequences \n")
        os.chdir(self.workdir)#Clean up dir moving
        pp = subprocess.call(["papara",
                              "-t", "random_resolve.tre",
                              "-s", "aln_ott.phy",
                              "-q", "{}.fasta".format(self.today),
                              "-n", papara_runname]) #FIx directory ugliness
        os.chdir('..')
        self.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.{}".format(self.workdir, papara_runname), schema="phylip")
        self.aln.taxon_namespace.is_mutable = False #This should enforce name matching throughout...
        self._query_seqs_aligned = 1
    def reconcile(self):
        """all missing data seqs are sneaking in, but from where?!"""
        if not self._query_seqs_aligned:
            self.align_query_seqs()
        aln_ids = set()
        for tax in self.aln:
            aln_ids.add(tax.label)
        assert aln_ids.issubset(self.otu_dict.keys())
        for key in  self.otu_dict.keys():
            if key not in aln_ids:
                sys.stderr.write("{} is in otu dict but not alignment\n".format(key))
        prune = []
        tmp_dict = {}
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        seq_len_cutoff = avg_seqlen*self.seq_len_perc
        for taxon, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(taxon.label)
            else:
                tmp_dict[taxon.label] = seq
                self.orig_seqlen.append(len(seq.symbols_as_string().translate(None, "-?")))
        self.aln = DnaCharacterMatrix.from_dict(tmp_dict)
        self.tre.prune_taxa_with_labels(prune)
        for tax in prune:
            del self.otu_dict[tax]
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("taxa pruned from tree and alignment due to excessive missing data\n")
            for tax in prune:
                fi.write("{}\n".format(tax))
            fi.close()
        self._reconciled = 1
    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new quesry seqs
        Just for placement, to use as starting tree."""
        sys.stdout.write("placing query sequences \n")
        os.chdir(self.workdir)
        p1 = subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                              "-f", "v",
                              "-s", "papara_alignment.extended",
                              "-t", "random_resolve.tre",
                              "-n", "PLACE"])
        placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                            schema="newick",
                            preserve_underscores=True)
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
        p2 = subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                              "-s", "papara_alignment.extended",
                              "-t", "place_resolve.tre",
                              "-p", "1",
                              "-n", "{}".format(self.today)])
        os.chdir('..')
        self._full_tree_est = 1
    def generate_streamed_alignment(self):
        """runs teh key steps and then replaces the tree and alignemnt with the expanded ones"""
        self.read_blast()
        if len(self.new_seqs) > 0:    
            self.remove_identical_seqs()
            self._write_files() #should happen before aligning in case of pruning
            self.write_query_seqs()
            self.align_query_seqs()
            self.reconcile()
            self.place_query_seqs()
            self.est_full_tree()
            self.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.workdir, self.today),
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
            self._write_files()
            if os.path.exists("{}/previous_run".format(self.workdir)):
                os.rename("{}/previous_run".format(self.workdir), "{}/previous_run{}".format(self.workdir, self.today))
            os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
            if os.path.exists("{}/last_completed_update".format(self.workdir)):
                os.rename(self.tmpfi,
                      "{}/last_completed_update".format(self.workdir))
            for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
            for filename in glob.glob('{}/papara*'.format(self.workdir)):
                os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
        else:
            sys.stdout.write("No new sequences found.")
        self.reset_markers()
        pickle.dump(self, open('{}/scrape.p'.format(self.workdir), 'wb'))
        pickle.dump(self.otu_dict, open('{}/otu_dict{}.p'.format(self.workdir, self.today), 'wb'))