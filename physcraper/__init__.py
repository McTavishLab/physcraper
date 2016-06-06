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


class ConfigObj(object):
    """Pulls out the configuration information from the config file and makes it easier to pass around and pickle."""
    def __init__(self, configfi):
        assert os.path.isfile(configfi)
        config = configparser.ConfigParser()
        config.read(configfi)
        self.e_value_thresh = config['blast']['e_value_thresh']
        self.hitlist_size = int(config['blast']['hitlist_size'])
        self.seq_len_perc = float(config['physcraper']['seq_len_perc'])
        self.get_ncbi_taxonomy = config['ncbi.taxonomy']['get_ncbi_taxonomy']
        self.ncbi_dmp = config['ncbi.taxonomy']['ncbi_dmp']
        self.phylesystem_loc = config['phylesystem']['location']
        self.ott_ncbi = config['ncbi.taxonomy']['ott_ncbi']




def generate_ATT_from_phylesystem(seqaln,
                                  mattype,
                                  workdir,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api'):
    """gathers together tree, alignment, and study info - forces names to otu_ids.
    Outputs AlignTreeTax object.
    an alignemnt, a
    Input can be either a study ID and tree ID from OpenTree"""
    #TODO CHECK ARGS
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
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
            print "{} doesn't have an otu id. What happened?".format(tax.label)#Forcing all spaces to underscore UGH
    otu_newick = tre.as_string(schema="newick")
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir) #newick should be bare, but alignement should be DNACharacterMatrix


def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            treefile,
                            otu_json,
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem.
    If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca."""
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH
    tre = Tree.get(path=treefile,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    with open(otu_json) as data_file:
        otu_dict = json.load(data_file)
    for tax in aln:
        assert tax.label in otu_dict
    tre = Tree.get(path=treefile,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    otu_newick = tre.as_string(schema="newick")
    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        ott_ids = [otu_dict[otu].get['^ot:ottId'] for otu in otu_dict]
        ott_mrca = get_mrca_ott(ott_ids)
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)


class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id """
    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir):
        self.aln = alignment
        self.tre = Tree.get(data=newick,
                            schema="newick",
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        self.otu_dict = otu_dict
        self.ps_otu = 1 #iterator for new otu IDs
        self._reconcile_names()
        self.workdir = workdir #TODO - is this where the workdir should live?
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca)
        self.ott_mrca = ingroup_mrca
        self.orig_seqlen = [1000, 1000, 1000] #FIXME
        self.gi_dict = {}
    def _reconcile_names(self):
        """This checks that the tree "original labels" from phylsystem
        align with those found in the taxonomy. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!"""
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
        missing = [i.label for i in self.aln.taxon_namespace if i.label not in treed_taxa]
        if missing:
            errmf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            raise ValueError(errm)
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
        return otu_id
    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """Papara is finicky about trees and needs phylip, this writes out needed files for papara
        (except query sequences)"""
        #CAN I even evaulte things in the function definitions?
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
    def write_otus(self, filename):
        """Writes out OTU dict as json"""
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)
    def write_labelled(self, label='^ot:ottTaxonName', treepath="labelled.tre", alnpath="labelled.fas"):
        """output tree and alignement with human readble labels
        Jumps through abunch of hoops to make labels unique.
        NOT MEMORY EFFICIENT AT ALL"""
        assert label in ['^ot:ottTaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
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
        tmp_tre.write(path="{}/{}".format(self.workdir, treepath),
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
        tmp_aln.write(path="{}/{}".format(self.workdir, alnpath),
                      schema="fasta")
#    def tidy_otu_dict():
#        #hmmm needs to strip out unused otus


#TODO... write as proper nexml?!


def prune_short(data_obj, min_seqlen=0):
    """Sometimes in the de-concatenating of the original alignment
    taxa with no sequence are generated.
    This gets rid of those from both the tre and the alignement. MUTATOR"""
    prune = []
    prune_orig = []
    tmp_dict = {}
    for taxon, seq in data_obj.aln.items():
        if len(seq.symbols_as_string().translate(None, "-?")) <= min_seqlen:
            prune.append(taxon.label)
            if data_obj.otu_dict[taxon.label].get(u'^ot:originalLabel'):
                prune_orig.append(data_obj.otu_dict[taxon.label].get(u'^ot:originalLabel'))
            else:
                prune_orig.append(taxon.label)
        else:
            tmp_dict[taxon.label] = seq
    data_obj.aln = DnaCharacterMatrix.from_dict(tmp_dict)
    data_obj.tre.prune_taxa_with_labels(prune)
    for tax in prune:
        del data_obj.otu_dict[tax]
    if prune:
        fi = open("{}/pruned_taxa".format(data_obj.workdir), 'w')
        fi.write("taxa pruned from tree and alignment due to excessive missing data\n")
        for tax in prune:
            fi.write("{}\n".format(tax))
        fi.close()


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
    try:
        mrca_node = tree_of_life.mrca(ott_ids=list(ott_ids), wrap_response=True)
    except  RuntimeError:
        sys.stderr.write("POST to get MRCA of ingroup failed - check internet connectivity, and or provide ingroup mrca OTT_ID, or check treemachine MRCA call\n")
        sys.exit()
    return mrca_node.nearest_taxon.ott_id


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
        fi.close()
        if os.path.isfile("{}/id_map.txt".format(workdir)): #todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]
    def add_gi(self, gi, tax_id):
        #add assert that it isn' already there?
        self.gi_ncbi_dict[gi] = tax_id
    def map_gi_ncbi(self, gi):
        """get the ncbi taxon id's for a gi input"""
        mapped_taxon_ids = open("{}/id_map.txt".format(self.workdir), "a")
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
            try:
                tax_id = int(subprocess.check_output(["bash", self.config.get_ncbi_taxonomy,
                                                      "{}".format(gi),
                                                      "{}".format(self.config.ncbi_dmp)]).split('\t')[1])
            except ValueError:
                sys.stderr.write("ncbi_dmp file needs to be updated. Do so and rerun.")
            mapped_taxon_ids.write("{}, {}\n".format(gi, tax_id))
            self.gi_ncbi_dict[gi] = tax_id
            assert tax_id  #if this doesn't work then the gi to taxon mapping needs to be updated - shouldhappen anyhow perhaps?!
        mapped_taxon_ids.close()
        return tax_id


#self.orig_seqlen = []

class PhyscraperScrape(object): #TODO do I wantto be able to instantiate this in a different way?!
    #set up needed variables as nones here?!
    #TODO better enforce ordering
    """This is the class that does the perpetual updating"""
    def __init__(self, data_obj, ids_obj, config_obj):
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
#        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "{}/{}.fasta".format(self.workdir, str(datetime.date.today()))
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
            last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
            today = str(datetime.date.today()).replace("-", "/")
            if abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days) > 14: #TODO make configurable
                equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi,
                                                               last_blast,
                                                               today)
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
                if not os.path.isfile(xml_fi):
                    try:
                        result_handle = NCBIWWW.qblast("blastn", "nt",
                                                       query,
                                                       entrez_query=equery,
                                                       hitlist_size=self.config.hitlist_size)
                        save_file = open(xml_fi, "w")
                        save_file.write(result_handle.read())
                        save_file.close()
                        self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        result_handle.close()
                    except (ValueError, URLError):
                        sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}".format(otu_id))
        self._blasted = 1
        return
    def read_blast(self, blast_dir=None):
        """reads in and prcesses the blast xml files"""
        if not blast_dir:
            blast_dir = self.blast_subdir
        if not self._blasted:
            self.run_blast()
        for taxon in self.data.aln:
            xml_fi = "{}/{}.xml".format(blast_dir, taxon.label)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                blast_records = NCBIXML.parse(result_handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if float(hsp.expect) < float(self.config.e_value_thresh):
                                if int(alignment.title.split('|')[1]) not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                    self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                    self.data.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
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
                    del self.data.otu_dict[label]
                    return
            else:
                if new_seq.find(inc_seq) != -1:#
                    if self.data.otu_dict[tax].get('^physcraper:status') == "original":
                        sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax))
                        seq_dict[label] = seq
                        return
                    else:
                        del seq_dict[tax]
                        self.data.tre.prune_taxa_with_labels(tax)
                        seq_dict[label] = seq
                        del self.data.otu_dict[tax]  #hmmmmmmmm?! CHECK ME!!
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
            log.write("{} new sequences added from genbank\n".format(len(self.new_seqs_otu_id)))
    def write_query_seqs(self):
        """writes out the query sequence file"""
        if not self._blast_read:
            self.read_blast()
        fi = open(self.newseqs_file, 'w')
        sys.stdout.write("writing out sequences\n")
        for otu_id in self.new_seqs_otu_id.keys():
            if otu_id not in self.data.aln: #new seqs only
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1
    def align_query_seqs(self, papara_runname="extended"):
        """runs papara on the tree, the alinment and the new query sequences"""
        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
                os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
        sys.stdout.write("aligning query sequences \n")
        self.data.write_papara_files()
        os.chdir(self.workdir)#Clean up dir moving
        pp = subprocess.call(["papara",
                              "-t", "random_resolve.tre",
                              "-s", "aln_ott.phy",
                              "-q", "{}.fasta".format(str(datetime.date.today())),
                              "-n", papara_runname]) #FIx directory ugliness
        os.chdir('..')
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = False #This should enforce name matching throughout...
        self._query_seqs_aligned = 1
    def reconcile(self):
        """all missing data seqs are sneaking in, but from where?!"""
        if not self._query_seqs_aligned:
            self.align_query_seqs()
        aln_ids = set()
        for tax in self.data.aln:
            aln_ids.add(tax.label)
        assert aln_ids.issubset(self.data.otu_dict.keys())
        for key in  self.data.otu_dict.keys():
            if key not in aln_ids:
                sys.stderr.write("{} is in otu dict but not alignment\n".format(key))
        prune = []
        tmp_dict = {}
        avg_seqlen = sum(self.data.orig_seqlen)/len(self.data.orig_seqlen)
        seq_len_cutoff = avg_seqlen*self.config.seq_len_perc
        for taxon, seq in self.data.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(taxon.label)
            else:
                tmp_dict[taxon.label] = seq
                self.data.orig_seqlen.append(len(seq.symbols_as_string().translate(None, "-?")))
        self.data.aln = DnaCharacterMatrix.from_dict(tmp_dict)
        self.data.tre.prune_taxa_with_labels(prune)
        for tax in prune:
            del self.data.otu_dict[tax]
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
                              "-n", "{}".format(str(datetime.date.today()))])
        os.chdir('..')
        self._full_tree_est = 1
    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignemnt with the expanded ones"""
        self.reset_markers()
        self.read_blast()
        if len(self.new_seqs) > 0:
            self.remove_identical_seqs()
            self.data.write_files() #should happen before aligning in case of pruning
            self.write_query_seqs()
            self.align_query_seqs()
            self.reconcile()
            self.place_query_seqs()
            self.est_full_tree()
            self.data.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.workdir, str(datetime.date.today())),
                                     schema="newick",
                                     preserve_underscores=True,
                                     taxon_namespace=self.data.aln.taxon_namespace)
            self.data.write_files()
            if os.path.exists("{}/previous_run".format(self.workdir)):
                os.rename("{}/previous_run".format(self.workdir), "{}/previous_run{}".format(self.workdir, str(datetime.date.today())))
            os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
            if os.path.exists("{}/last_completed_update".format(self.workdir)):
                os.rename(self.tmpfi, "{}/last_completed_update".format(self.workdir))
            for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
            for filename in glob.glob('{}/papara*'.format(self.workdir)):
                os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
        else:
            sys.stdout.write("No new sequences found.")
        self.reset_markers()
        pickle.dump(self, open('{}/scrape.p'.format(self.workdir), 'wb'))
        pickle.dump(self.data.otu_dict, open('{}/otu_dict{}.p'.format(self.workdir, str(datetime.date.today())), 'wb'))

'''
class PhyscraperSetup(object):
    """This needs to vet the inputs, standardize names, and prep for
    first run through of blasting. Does not blast itself!"""
    _found_mrca = 0
    _id_dicts = 0
    _study_get = 0
    _phylesystem = 0
    def __init__(self,
                 AlignTreeTax,
                 IdDicts,
                 study_info,
                 runname,
                 configfi='/home/ejmctavish/projects/otapi/physcraper/config'):
        """initialized object, most attributes generated through self._checkArgs using config file."""
#        self.study_id = study_info.study_id
 #       self.tree_id = study_info.tree_id
 #       self.seqaln = study_info.seqaln
 #       self.mattype = study_info.mattype
 #       self.treefile = study_info.treefile
 #       self.otu_dict = study_info.otu_dict
        self.runname = runname
#        self._get_study()
#        self.ott_to_ncbi = {}
#        self.ott_to_name = {}
#        self.ncbi_to_ott = {}
#        self.phy = None
#        self.gi_ncbi_dict = {}
        self.mrca_ott = None
#        self.aln = None
#        self.orig_aln = None
#        self.tre = None
        self.mrca = get_mrca_ott(self.ott_ids)
        mrca_ncbi = self.ott_to_ncbi[self.mrca_ott]
#        self.orig_seqlen = None
#        try:
#            self._reconcile_names()
#        except:
#            print "name reconcillation error"
#        self.workdir = runname
        self.data = AlignTreeTax
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
#    def write_files(self):
#        self.tre.write(path="{}/original_pruned.tre".format(self.workdir),
#                       schema="newick",
#                       unquoted_underscores=True)
#        self.aln.write(path="{}/original_pruned.fas".format(self.workdir),
#                       schema="fasta")
#    def __getstate__(self):
#        """hacky way to fix pickling"""
#        result = self.__dict__.copy()
#        del result['config']
#        del result['phy']
#        return result
#    def setup_physcraper(self):
#        """run the key steps"""
#        self._read_config()
#        self._get_mrca()
#        self._reconcile_names()
#        self._prune()
#        self.write_files()
# I think it should write out to files so that the PHyscarpe scrape object can just read the aln etc from those, maybe the shaould be where the OTU dict comes from too...
#Hmmm.


#    def _write_files(self):
#        """Outputs both the streaming files and a ditechecked"""
#        #First write rich annotation json file with everything needed for later?
#        self.tre.resolve_polytomies()
#        self.tre.deroot()
#        tmptre = self.tre.as_string(schema="newick",
#                                    unquoted_underscores=True,
#                                    suppress_rooting=True)
#        tmptre = tmptre.replace(":0.0;",";")#Papara is diffffffficult about root
#        fi = open("{}/random_resolve.tre".format(self.workdir), "w")
#        fi.write(tmptre)
#        fi.close()
#        self.aln.write(path="{}/aln_ott.phy".format(self.workdir),
#                       schema="phylip")
#        self.aln.write(path="{}/aln_ott.fas".format(self.workdir),
#                       schema="fasta")
#        self.tre.write(path="{}/random_resolve{}.tre".format(self.workdir, self.today),
#                       schema="newick",
#                       unquoted_underscores=True)
#        self.aln.write(path="{}/aln_ott{}.fas".format(self.workdir, self.today),
#                      schema="fasta")

#        if setup_obj:#this should need to happen ony the first time the scrape object is run.
#            self.runname = 
#            
#            self.e_value_thresh = deepcopy(setup_obj.e_value_thresh)
#            self.hitlist_size = deepcopy(setup_obj.hitlist_size)
#            self.seq_len_perc = deepcopy(setup_obj.seq_len_perc)
#            self.get_ncbi_taxonomy = deepcopy(setup_obj.get_ncbi_taxonomy)
#            self.ncbi_dmp = deepcopy(setup_obj.ncbi_dmp)
#            self.orig_seqlen = deepcopy(setup_obj.orig_seqlen)
        self.workdir = data_obj.workdir #TODO potnetial from clash between ids dict data obj
'''