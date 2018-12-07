#!/usr/bin/env python
import os
import sys
import subprocess
import csv
import pickle
import random

from copy import deepcopy
from dendropy import Tree, DnaCharacterMatrix
from Bio import Entrez

import physcraper
# from __init__ import cd

if sys.version_info < (3, ):
    from urllib2 import HTTPError
else:
    from urllib.error import HTTPError

"""Code used to concatenate different single PhyScraper runs into a concatenated one.
"""

# def remove_leaf(tre, leaf):
#     """ Removes a taxon from tre.
#     Does not work the way it is intended, as the flexibility of leaf gets lost.
#     """
#     # not used somewhere
#     tre.prune_taxa([leaf])
#     tre.prune_taxa_with_labels([leaf.taxon])
#     tre.prune_taxa_with_labels([leaf])


def remove_aln_tre_leaf(scrape):
    """Attempt to remove all occurrences in aln and tre of otu,
     that were removed sometime in the single runs, but kind of sneak in.
    """
    aln_ids = set()
    for tax in scrape.data.aln:
        aln_ids.add(tax.label)
    assert aln_ids.issubset(scrape.data.otu_dict.keys())

    treed_taxa = set()
    for leaf in scrape.data.tre.leaf_nodes():
        treed_taxa.add(leaf.taxon)
    for leaf in scrape.data.tre.leaf_nodes():
        if leaf.taxon not in aln_ids:
            scrape.data.tre.prune_taxa([leaf])
            scrape.data.tre.prune_taxa_with_labels([leaf.taxon])
            scrape.data.tre.prune_taxa_with_labels([leaf])
            treed_taxa.remove(leaf.taxon)
    assert treed_taxa.issubset(aln_ids)
    return scrape


def add_to_del_acc(del_acc, gene, spn, random_gen):
    """Adds gi number to del_acc.
    Del_acc is used to remove gi's from tmp_dict, so that they will
    not be added to the concat dict twice.
    """
    spn_ = spn.replace(" ", "_")
    if gene in del_acc.keys():
        if spn_ not in del_acc[gene].keys():
            del_acc[gene][spn] = random_gen
        else:
            del_acc[gene][spn_.format("_", " ")] = random_gen
    else:
        del_acc[gene] = {spn: random_gen}
    return del_acc


class Concat(object):
    """Combines several physcraper runs into a concatenated alignment and calculates a phylogeny.
    There are two options available, either data will be concatenated by random (per taxon name) or the 
    user provides a file which say, which sequences shall be concatenated.

    User need to make sure, that there are at least some overlapping lineages.
    Do not concatenate data from the same loci (if you want to expand an alignment, run physcraper!).

    To build the class the following is needed:
        workdir_comb: the path to your directory where the data shall be stored
        email: your email address, currently used to retrieve missing taxon information

    During the initializing process the following self objects are generated:
        self.workdir: the path to your directory
        self.sp_acc_comb: dictonary - similar to otu_dcit of Physcraper class
            key: taxon name
            value: dictionary:
                key: gene name
                value: dictionary:
                    key: unique identifier of concat_class
                    value: dictionary - key-value-pairs:
                        "unique_id": unique id - either user name or genbank accession number
                        "seq": corresponding sequence
                        "spn": taxon name
                        "original_PS_id": otu_id from single-gene run
                        "concat:status": "single run"/ "concatenated"
                        "new tipname": taxon name plus number if there are more than a single
                                        concatenated sequence per taxon
        self.single_runs: dictonary
            key: gene name, as provided by user in input
            value: file containing the single gene Physcraper run, loaded from pickle
        self.sp_counter: dictonary
            key: taxon name
            value: dictionary
                key: gene name
                value: number of sequences present for the given gene and taxon
        self.email: email
        self.comb_seq: dictonary
            key: gene name
            value: dictionary:
                key: taxon name
                value: sequence
        self.comb_acc: dictonary # !!! Note can be easily combined with comb_seq
            key: gene name
            value: dictionary:
                key: taxon name
                value: unique identifier of concat_class
        self.aln_all: dictonary
            key: numbers from 1 to amount of loci concatenated
            value: dendropy aln including empty seq
        self.num_of_genes = number corresponding to the number of genes that shall be concatenated
        self.genes_present = list of the gene names that shall be concatenated
        self.tre_as_start = phylogeny used as starting tree for concatenated run, is the one with most tips present
        self.tre_start_gene = corresponding name to tre_as_start
        self.short_concat_seq = list of taxa that have only few sequences information/short genes
        self.concat_tips: dictonary
            key: otu.label from self.tre_as_start.taxon_namespace
            value: taxon name
        self.concatfile = path to file if user supplied concatenation file is used
        self.concatenated_aln = concatenated alignment
        self.tmp_dict = subset of self.sp_acc_comb
        self.part_len = holds sequence partition position to write the partitioning file # ! TODO MK: might not need to be self
    """

    def __init__(self, workdir_comb, email):
        # super(PhyscraperScrape, self).__init__()
        self.workdir = workdir_comb
        self.sp_acc_comb = {}
        self.single_runs = {}
        self.sp_counter = {}
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.email = email
        self.comb_seq = {}
        self.comb_acc = {}
        self.aln_all = {}
        self.num_of_genes = 0
        self.genes_present = []
        self.tre_as_start = None
        self.tre_start_gene = None
        self.short_concat_seq = None
        self.concatfile = None
        self.part_len = None
        self.concatenated_aln = None
        self.tmp_dict = None
        self.concat_tips = {}

    def load_single_genes(self, workdir, pickle_fn, genename):
        """Load PhyScraper class objects and make a single dict per run.

        Removes abandoned nodes first.

        :param workdir: directory of single gene run
        :param pickle_fn: path to pickled file of the Physcraper run
        :param genename: string, name for locus provided by user
        :return: self.single_runs
        """
        physcraper.debug("load_single_genes: {}".format(genename))
        scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn), "rb"))
        scrape = remove_aln_tre_leaf(scrape)
        self.single_runs[genename] = deepcopy(scrape)
        return

    def combine(self):
        """Combines several PhyScraper objects to make a concatenated run dict.

        Is a wrapper function around make_concat_id_dict(). It produces the parameters needed for the function.
        """
        physcraper.debug("combine")
        self.num_of_genes = len(self.single_runs)
        concat_id_counter = 1
        for genename in self.single_runs:
            self.genes_present.append(genename)
            for otu in self.single_runs[genename].data.aln.taxon_namespace:
                concat_id = "concat_{}".format(concat_id_counter)
                self.make_concat_id_dict(otu.label, genename, concat_id)
                concat_id_counter += 1
        return

    def make_concat_id_dict(self, otu, genename, concat_id):
        """Makes a concat_id entry with all information

        Note: has test

        :param otu: otu_id
        :param genename: name of single gene run
        :param concat_id: unique identifier in the concat class
        :return: modified self.sp_acc_comb
        """
        data = self.single_runs[genename].data.otu_dict[otu]
        seq = str(self.single_runs[genename].data.aln[otu])
        spn = None
        if "^ot:ottTaxonName" in data:
            spn = self.get_taxon_info("^ot:ottTaxonName", data)
            if spn not in self.sp_acc_comb:
                self.sp_acc_comb[spn] = {}
            if genename not in self.sp_acc_comb[spn]:
                self.sp_acc_comb[spn][genename] = {}
        elif "^user:TaxonName" in data:
            spn = self.get_taxon_info("^user:TaxonName", data)
            if spn not in self.sp_acc_comb:
                self.sp_acc_comb[spn] = {}
            if genename not in self.sp_acc_comb[spn]:
                self.sp_acc_comb[spn][genename] = {}
        else:
            # we should never get here....
            physcraper.debug("THERE IS A SERIOUS PROBLEM....")
        assert spn is not None
        if concat_id not in self.sp_acc_comb[spn][genename]:
            if "^ncbi:accession" in data:
                unique_id = data["^ncbi:accession"]
            elif u"^ot:originalLabel" in data:
                unique_id = data[u"^ot:originalLabel"]
            concat_dict = {
                "unique_id": unique_id,
                "seq": seq,
                "spn": spn,
                "original_PS_id": otu,
                "concat:status": "single run",
            }
            self.sp_acc_comb[spn][genename][concat_id] = concat_dict
        else:
            physcraper.debug(
                "something goes wrong, you should not try to add the same id several times...."
            )
        if concat_dict["spn"] is None:
            # we should never get here....
            sys.stderr.write(
                "There is no species name for the seq. Do not know how to concatenate then. "
                "Please remove seq from aln: {}.".format(data["^ncbi:accession"])
            )
            physcraper.debug("THERE IS A SERIOUS PROBLEM....spn is none")
            spn = self.get_taxon_info("^ot:ottTaxonName", data)
            self.sp_acc_comb[spn] = self.sp_acc_comb[unique_id]
            del self.sp_acc_comb[unique_id]

    def get_taxon_info(self, key, data):
        """If there are no taxon information (for what ever reason) try again to obtain sp names.

        If the key is not part of data, it will get the name through a web query using the GI number.

        :param key: key of otu_dict/data that shall contain the taxon name, e.g.^ot:ottTaxonName
        :param data: otu_dict entry from single gene physcraper run
        :return: taxon name
        """
        # debug("get_rank_info")
        tax_name = None
        if key in data:
            if data[key] is None:
                if "^ncbi:accession" in data:
                    gb_id = data["^ncbi:accession"]
                Entrez.email = self.email
                tries = 5
                for i in range(tries):
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=gb_id, retmode="xml")
                    except (IndexError, HTTPError) as err:
                        sys.stderr.write(err)
                        if i < tries - 1:  # i is zero indexed
                            continue
                        else:
                            raise
                    break
                read_handle = Entrez.read(handle)[0]
                tax_name = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value']
            else:
                tax_name = data[key]
        assert tax_name is not None
        tax_name = tax_name.replace("_", " ")
        tax_name = tax_name.replace(".", "").replace("'", "")
        tax_name = tax_name.encode("ascii")
        return tax_name

    def sp_seq_counter(self):
        """Counts how many seq per sp and genes there are -is used by sp_to_keep.

        Note: has test

        :return: builds self.sp_counter
        """
        physcraper.debug("sp_seq_counter")
        for spn in self.sp_acc_comb:
            tmp_gene = deepcopy(self.genes_present)
            for gene in self.sp_acc_comb[spn]:
                tmp_gene.remove(gene)
                spn_new = spn.replace(" ", "_")
                if spn_new in self.sp_counter:
                    self.sp_counter[spn_new][gene] = len(self.sp_acc_comb[spn][gene])
                else:
                    self.sp_counter[spn_new] = {gene: len(self.sp_acc_comb[spn][gene])}
            for item in tmp_gene:
                if spn_new in self.sp_counter:
                    self.sp_counter[spn_new][item] = 0
                else:
                    self.sp_counter[spn_new] = {item: 0}
        physcraper.debug(self.sp_counter)

    def get_largest_tre(self):
        """Find the single gene tree with the most tips, which will be used as
        starting tree for concat phylo reconstruction.
        """
        physcraper.debug("get_largest_tre")
        first = True
        len_all_taxa = {}
        for gene in self.single_runs:
            len_aln_taxa = len(self.single_runs[gene].data.aln.taxon_namespace)
            len_all_taxa[gene] = len_aln_taxa
        len_max = 0
        gene_max = 0
        for gene, len_item in len_all_taxa.items():
            if first:
                len_max = len_item
                gene_max = gene
                assert len_max != 0
                assert gene_max != 0
                first = False
            if len_item > len_max:
                len_max = len_item
                gene_max = gene
        self.tre_as_start = self.single_runs[gene_max].data.tre
        self.tre_start_gene = gene_max

    def make_sp_gene_dict(self):
        """Is the build around to make the dicts that are used to make it into a dendropy aln
        """
        physcraper.debug("make_sp_gene_dict")
        if self.concatfile is not None:
            self.user_defined_concat()
        else:
            sp_to_keep = self.sp_to_keep()
            self.tmp_dict = deepcopy(self.sp_acc_comb)
            while len(self.tmp_dict.keys()) >= 1:
                del_acc = {}
                for spn in self.tmp_dict.keys():
                    sp_to_keep_list = sp_to_keep.keys()
                    if spn.replace(" ", "_") in sp_to_keep_list:
                        tmp_gene = deepcopy(self.genes_present)
                        for gene in self.tmp_dict[spn]:
                            tmp_gene.remove(gene)
                            del_acc = self.select_rnd_seq(spn, gene, del_acc)
                        for item in tmp_gene:
                            self.make_empty_seq(spn, item)
                        self.rm_rnd_sp(del_acc)
                        del self.tmp_dict[spn]
                    else:
                        for gene in self.tmp_dict[spn]:
                            del_acc = self.select_rnd_seq(spn, gene, del_acc)
                        self.rm_rnd_sp(del_acc)
                    self.rm_empty_spn_entries(del_acc)
        self.rename_drop_tips()

    def rename_drop_tips(self):
        """ Removes tips from tre as start that are not present in the concatenated aln
        and renames tips that are present.
        """
        physcraper.debug("rename_drop_tips")
        # leaf.taxon is never in concat_tips
        for leaf in self.tre_as_start.leaf_nodes():
            if leaf.taxon.label not in self.concat_tips.keys():
                self.tre_as_start.prune_taxa([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.label])
                self.tre_as_start.prune_taxa_with_labels([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
            else:
                for otu in self.concat_tips.keys():
                    if otu == leaf.taxon.label:
                        leaf.taxon.label = self.concat_tips[otu]

    def sp_to_keep(self):
        """Uses the sp_counter to make a list of sp that should be kept in concatenated alignment,
        because they are the only representative of the sp.

        Note: has test

        :return: dictionary with taxon name and number saying how many genes are missing
        """
        physcraper.debug("sp to keep")
        sp_to_keep = {}
        for spn in self.sp_counter:
            seq_counter = True
            not_present = 0
            for gene in self.sp_counter[spn]:
                if self.sp_counter[spn][gene] == 0:
                    seq_counter = False
                    not_present += 1
            if not seq_counter:
                sp_to_keep[spn] = not_present
        # physcraper.debug(sp_to_keep)
        return sp_to_keep

    def select_rnd_seq(self, spn, gene, del_acc):
        """Select a random seq from spn and gene to combine it with a random other one from another gene,
        but same spn. Is used if the user does not give a concatenation input file.

        Note: has test

        :param spn: taxon name
        :param gene:  gene name
        :param del_acc: dictionary that contains gene name: dict(spn: concat_id of random seq)
        :return: del_acc
        """
        physcraper.debug("select_rnd_seq")
        count = 2
        random_gen = random.choice(list(self.tmp_dict[spn][gene]))
        self.sp_acc_comb[spn][gene][random_gen]["concat:status"] = "used in concat"
        seq = str(self.tmp_dict[spn][gene][random_gen]["seq"])
        spn_ = spn.replace(" ", "_")
        spn_ = spn_.replace(".", "").replace("'", "")
        if gene in self.comb_seq.keys():
            if spn_ not in self.comb_seq[gene].keys():
                self.comb_seq[gene][spn_] = seq
                if gene in self.comb_acc:
                    self.comb_acc[gene][spn_] = random_gen
                else:
                    self.comb_acc[gene] = {spn_: random_gen}
                if gene in del_acc.keys():
                    if spn_ not in del_acc[gene].keys():
                        del_acc[gene][spn] = random_gen
                else:
                    del_acc[gene] = {spn: random_gen}
            else:
                spn_new = "{}_{}".format(spn_, count)
                while spn_new in self.comb_seq[gene].keys():
                    count += 1
                    spn_new = "{}_{}".format(spn_, count)
                self.comb_seq[gene][spn_new] = seq
                self.comb_acc[gene][spn_new] = random_gen
                self.sp_acc_comb[spn][gene][random_gen]["new tipname"] = spn_new
                if gene in del_acc.keys():
                    if spn_ not in del_acc[gene].keys():
                        del_acc[gene][spn] = random_gen
                    else:
                        del_acc[gene][spn] = random_gen
                else:
                    del_acc[gene] = {spn: random_gen}
        else:
            self.comb_seq[gene] = {spn_: seq}
            self.comb_acc[gene] = {spn_: random_gen}
            if gene in del_acc.keys():
                if spn_ not in del_acc[gene].keys():
                    del_acc[gene][spn] = random_gen
                else:
                    del_acc[gene] = {spn: random_gen}
            else:
                del_acc[gene] = {spn: random_gen}
        self.otu_to_spn(spn, gene, del_acc[gene][spn])
        return del_acc

    def otu_to_spn(self, spn_, gene, random_gen):
        """ Makes a dict that contains the original tip labels from the starting tree,
        and the name it needs to have (no name duplicates allowed for many phylogenetic programs).
        This dict will be used in rename_drop_tips to rename or remove the tips.
        As such it is a helper function to produce the correct starting tree, for the concatenation run.

        :param spn_: species name for concat
        :param gene: gene
        :param random_gen: the corresponding otu
        :return: self.concat_tips
        """
        if self.tre_start_gene == gene:
            spn = spn_.replace("_", " ")
            former_otu = self.sp_acc_comb[spn][gene][random_gen]["original_PS_id"]
            for otu in self.tre_as_start.taxon_namespace:
                if otu.label == former_otu:
                    if "new tipname" in self.sp_acc_comb[spn][gene][random_gen]:
                        spn_ = self.sp_acc_comb[spn][gene][random_gen]["new tipname"]
                    self.concat_tips[otu.label] = spn_
        return self.concat_tips

    def make_empty_seq(self, spn, gene):
        """Is called when there is no seq for one of the genes, but otu shall be kept in aln.
        Dendropy needs same taxon_namespace and number otu's for concatenation. It will just make an empty sequence of
        the same length.
        """
        # physcraper.debug("make_empty_seq")
        len_gene_aln = 0
        for tax, seq in self.single_runs[gene].data.aln.items():
            len_gene_aln = len(seq)
            break
        assert len_gene_aln != 0
        empty_seq = "?" * len_gene_aln
        if gene in self.comb_seq:
            self.comb_seq[gene][spn.replace(" ", "_")] = empty_seq
        else:
            self.comb_seq[gene] = {spn.replace(" ", "_"): empty_seq}

    def rm_rnd_sp(self, del_acc):
        """Removes the random selected seq from the tmp_dict, so that it cannot be selected again.
        """
        physcraper.debug("rm_rnd sp")
        for spn2 in self.tmp_dict:
            for gene2 in self.tmp_dict[spn2]:
                if gene2 in del_acc:
                    if spn2 in del_acc[gene2]:
                        key = del_acc[gene2][spn2]
                        if key in self.tmp_dict[spn2][gene2]:
                            del self.tmp_dict[spn2][gene2][key]

    def rm_empty_spn_entries(self, del_acc):
        """Removes keys from tmp dict, if the key/sp has no value anymore. Helper function.
        """
        physcraper.debug("rm_empty_spn_entries")
        del_sp = None
        for spn2 in self.tmp_dict:
            for gene2 in self.tmp_dict[spn2]:
                if gene2 in del_acc:
                    if spn2 in del_acc[gene2]:
                        if len(self.tmp_dict[spn2][gene2]) == 0:
                            del_sp = spn2
        if del_sp is not None:
            for item in self.sp_acc_comb[del_sp]:
                for otu in self.sp_acc_comb[del_sp][item]:
                    if self.sp_acc_comb[del_sp][item][otu]["concat:status"] != "used in concat":
                        self.sp_acc_comb[del_sp][item][otu][
                            "concat:status"] = "deleted, because not enough seq are present"
            del self.tmp_dict[del_sp]

    def make_alns_dict(self):
        """Makes dendropy aln out of dict self.comb_seq for all genes.
        """
        physcraper.debug("make_alns_dict")
        firstelement = True
        count = 0
        for gene in self.comb_seq.keys():
            if count == 0:
                len1 = len(self.comb_seq[gene].keys())
                len2 = len1
                count = 1
            else:
                len2 = len(self.comb_seq[gene].keys())
            assert len1 == len2
        for gene in self.comb_seq.keys():
            if firstelement:
                aln1 = DnaCharacterMatrix.from_dict(self.comb_seq[gene])
                firstelement = False
                self.aln_all[count] = aln1
                # aln1.write(path="{}/aln_0.fas".format(self.workdir),
                           # schema="fasta")
            else:
                aln = DnaCharacterMatrix.from_dict(self.comb_seq[gene], taxon_namespace=aln1.taxon_namespace)
                self.aln_all[count] = aln
                # aln.write(path="{}/aln_{}.fas".format(self.workdir, count),
                          # schema="fasta")
            count += 1

    def concatenate_alns(self):
        """Concatenate all alns into one aln.
        """
        physcraper.debug("concat alns")
        count = 0
        for gene in self.aln_all:
            if count == 0:
                aln1 = self.aln_all[gene]
                aln1.write(path="{}/aln1.fas".format(self.workdir), schema="fasta")
                count = 1
            else:
                aln2 = self.aln_all[gene]
                count += 1
                aln2.write(path="{}/aln{}.fas".format(self.workdir, count), schema="fasta")
                assert aln1.taxon_namespace == aln2.taxon_namespace
                aln1 = DnaCharacterMatrix.concatenate([aln1, aln2])
        aln1.write(path="{}/concat.fas".format(self.workdir),
                   schema="fasta")
        self.concatenated_aln = aln1

    def get_short_seq_from_concat(self, percentage=0.37):
        """Finds short sequences, all below a certain threshold will be removed,
        to avoid having really low coverage in the aln. Default = 0.37.

        Note percentage is a bit misleading, the cutoff is 37% of the whole concatenated
        alignment, but the sequences length is calculated without gaps present.
        The default is so low, as I want to keep taxa that have only a single locus
        and which is not the longest among the loci within the aln.
        """
        physcraper.debug("get_short_seq_from_concat")
        seq_len = {}
        num_tax = 0
        for tax, seq in self.concatenated_aln.items():
            seq = seq.symbols_as_string().replace("-", "").replace("?", "")
            seq_len[tax] = len(seq)
            num_tax += 1
        total_len = 0
        for tax, seq in self.concatenated_aln.items():
            total_len = len(seq)
            break
        assert total_len != 0
        min_len = total_len * percentage
        prune_shortest = []
        for tax, len_seq in seq_len.items():
            if len_seq < min_len:
                prune_shortest.append(tax)
        self.short_concat_seq = prune_shortest

    def remove_short_seq(self):
        """Removes short seq that were found with get_short_seq
        and write it to file.
        """
        physcraper.debug("remove_short_seq")
        self.concatenated_aln.remove_sequences(self.short_concat_seq)
        for leaf in self.tre_as_start.leaf_nodes():
            for tax in self.short_concat_seq:
                if tax.label == leaf.taxon.label.replace(" ", "_"):
                    self.tre_as_start.prune_taxa([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.label])
                    self.tre_as_start.prune_taxa_with_labels([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                    self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
                else:
                    leaf.taxon.label = leaf.taxon.label.replace(" ", "_")
        tre_as_start_str = self.tre_as_start.as_string(schema="newick",
                                                       # preserve_underscores=True,
                                                       unquoted_underscores=True,
                                                       suppress_rooting=True)
        fi = open("{}/{}".format(self.workdir, "starting_red.tre"), "w")
        fi.write(tre_as_start_str)
        fi.close()
        for tax in self.concatenated_aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")
        self.concatenated_aln.write(path="{}/{}".format(self.workdir, "concat_red.fasta"), schema="fasta")
        tre_ids = set()
        for tax in self.tre_as_start.taxon_namespace:
            tre_ids.add(tax.label)
        aln_ids = set()
        for tax in self.concatenated_aln.taxon_namespace:
            aln_ids.add(tax.label)

    def make_concat_table(self):
        """ Makes a table that shows which gi numbers were combined.
        """
        genel = []
        spn_l = {}
        for gene in self.comb_acc:
            genel.append(gene)
            for spn in self.comb_acc[gene]:
                concat_id = self.comb_acc[gene][spn]
                multiname = False
                if spn.split("_")[-1].isdigit():
                    tmp_spn = spn[:-2]
                    tmp_spn = tmp_spn.replace("_", " ")
                    multiname = True
                spn = spn.replace("_", " ")
                if multiname:
                    if tmp_spn in self.sp_acc_comb.keys():
                        unique_id = self.sp_acc_comb[tmp_spn][gene][concat_id]["unique_id"]
                    else:
                        unique_id = self.sp_acc_comb[spn][gene][concat_id]["unique_id"]
                else:
                    unique_id = self.sp_acc_comb[spn][gene][concat_id]["unique_id"]
                if spn in spn_l.keys():
                    spn_l[spn].append(unique_id)
                else:
                    spn_l[spn] = [unique_id]
        with open("{}/concatenation.csv".format(self.workdir), "w") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(genel)
            for key, value in spn_l.items():
                writer.writerow([key, value])

    def write_partition(self):
        """Write the partitioning file for RAxML.
        """
        physcraper.debug("write_partition")
        count = 0
        len_gene = 0
        for gene in self.single_runs:
            for tax, seq in self.single_runs[gene].data.aln.items():
                len_gene = len(seq.symbols_as_string())
                break
            if count == 0:
                with open("{}/partition".format(self.workdir), "w") as partition:
                    partition.write("DNA, {} = 1-{}\n".format(gene, len_gene))
                self.part_len = len_gene
                count = 1
            else:
                start = self.part_len + 1
                end = self.part_len + len_gene
                self.part_len = self.part_len + len_gene
                with open("{}/partition".format(self.workdir), "a") as partition:
                    partition.write("DNA, {} = {}-{}\n".format(gene, start, end))

    def place_new_seqs(self):
        """Places the new seqs (that are only found in loci which is not the starting tree)
        onto one of the single run trees.
        """
        physcraper.debug("place_new_seqs")
        if len(self.concatenated_aln.taxon_namespace)-len(self.short_concat_seq) > len(self.tre_as_start.leaf_nodes()):
            if os.path.exists("RAxML_labelledTree.PLACE"):
                os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")

            with physcraper.cd(self.workdir):
                # cwd = os.getcwd()
                # os.chdir(self.workdir)

                physcraper.debug("make place-tree")
                try:
                    num_threads = int(self.config.num_threads)
                    print(num_threads)
                    subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                     "-f", "v", "-q", "partition",
                                     "-s", "concat_red.fasta",
                                     "-t", "starting_red.tre",
                                     "-n", "PLACE"])
                except:
                    subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                     "-f", "v", "-q", "partition",
                                     "-s", "concat_red.fasta",
                                     "-t", "starting_red.tre",
                                     "-n", "PLACE"])
                # os.chdir(cwd)
            physcraper.debug("read place tree")
            placetre = Tree.get(path="{}/starting_red.tre".format(self.workdir),
                                schema="newick",
                                preserve_underscores=True,
                                suppress_internal_node_taxa=True, suppress_leaf_node_taxa=True)
            physcraper.debug("resolve polytomies")
            placetre.resolve_polytomies()
            placetre.write(path="{}/place_resolve.tre".format(self.workdir), schema="newick", unquoted_underscores=True)

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree.
        """
        physcraper.debug("run full tree")
        with physcraper.cd(self.workdir):
            # cwd = os.getcwd()
            # os.chdir(self.workdir)
            if os.path.exists("place_resolve.tre"):
                starting_fn = "place_resolve.tre"
            else:
                starting_fn = "starting_red.tre"
            if os.path.exists("concat_red.fasta.reduced"):
                aln = "concat_red.fasta.reduced"
                partition = "partition.reduced"
            else:
                aln = "concat_red.fasta"
                partition = "partition"
            try:
                num_threads = int(self.config.num_threads)
                print(num_threads)
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                             "-s", aln, "--print-identical-sequences",
                             "-t", "{}".format(starting_fn),
                             "-p", "1", "-q", partition,
                             "-n", "concat"])
            except:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-s", aln, "--print-identical-sequences",
                             "-t", "{}".format(starting_fn),
                             "-p", "1", "-q", partition,
                             "-n", "concat"])
            # os.chdir(cwd)

    def calculate_bootstrap(self):
        """Calculate bootstrap and consensus trees.
        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        """
        with physcraper.cd(self.workdir):
            # os.chdir(self.workdir)
            if os.path.exists("concat_red.fasta.reduced"):
                aln = "concat_red.fasta.reduced"
                partition = "partition.reduced"
            else:
                aln = "concat_red.fasta"
                partition = "partition"
            # run bootstrap
            # make bipartition tree
            # is the -f b command
            # -z specifies file with multiple trees
            try:
                num_threads = int(self.config.num_threads)
                print(num_threads)
                # subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                #                  "-s", aln,  "-q", partition,
                #                  # "-t", "place_resolve.tre",
                #                  "-p", "1", "-b", "1", "-#", "autoMRE",
                #                  "-n", "autoMRE"])
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", aln, "-q", partition,
                                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                                 "-n", "autoMRE_fa"])
                # strict consensus:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-J", "STRICT",
                                 "-z", "RAxML_bootstrap.autoMRE_fa",
                                 "-n", "StrictCon"])
                # majority rule:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-J", "MR",
                                 "-z", "RAxML_bootstrap.autoMRE_fa",
                                 "-n", "MR"])
            except:
                # subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                #                  "-s", aln,  "-q", partition,
                #                  # "-t", "place_resolve.tre",
                #                  "-p", "1", "-b", "1", "-#", "autoMRE",
                #                  "-n", "autoMRE"])
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", aln, "-q", partition,
                                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                                 "-n", "autoMRE_fa"])
                # strict consensus:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-J", "STRICT",
                                 "-z", "RAxML_bootstrap.autoMRE_fa",
                                 "-n", "StrictCon"])
                # majority rule:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-J", "MR",
                                 "-z", "RAxML_bootstrap.autoMRE_fa",
                                 "-n", "MR"])

    def user_defined_concat(self):
        """If a user gave an input file to concatenate data. Fills in the data for self.comb_seq, self.comb_acc
        (is the replacement function for select_rnd_seq).
        """
        physcraper.debug("user_defined_concat")
        with open("{}/{}".format(self.workdir, self.concatfile), mode="r") as infile:
            reader = csv.reader(infile)
            sp_concat = dict((rows[0], rows[1]) for rows in reader)
        for otu in sp_concat.keys():
            global_spn = None
            concat_l = sp_concat[otu]
            if concat_l[:1] == "[":
                concat_l = concat_l[1:-1]
            concat_l = concat_l.split(", ")
            for item in concat_l:
                gene_l = []
                if item[:1] == "'":
                    item = item[1:-1]
                item = item.encode("utf-8")
                for gene in self.single_runs:
                    spn = None
                    for key, val in self.single_runs[gene].data.otu_dict.items():
                        if item.isdigit():
                            if "^ncbi:accession" in val:
                                if item == val["^ncbi:accession"]:
                                    spn = val["^ot:ottTaxonName"]
                                    gene_l.append(gene)
                        else:
                            if "^ncbi:accession" in val:
                                if item == val["^ncbi:accession"]:
                                    spn = val["^ot:ottTaxonName"]
                                    gene_l.append(gene)
                            elif u"^ot:originalLabel" in val:
                                if item == val[u"^ot:originalLabel"]:
                                    spn = val["^ot:ottTaxonName"]
                                    gene_l.append(gene)
                        if spn is not None:
                            global_spn = spn.replace(".", "").replace("'", "")
                            spn = spn.replace(".", "").replace("'", "")
                            for key2, val2 in self.sp_acc_comb[spn][gene].items():
                                cond = False
                                if len(item.split(".")) >= 2 and val2["unique_id"] == item:
                                    cond = True
                                else:
                                    if val2["unique_id"] == item:
                                        cond = True
                                if cond:
                                    concat_id = key2
                                    self.sp_acc_comb[spn][gene][concat_id]["concat:status"] = "used in concat"
                                    seq = str(self.sp_acc_comb[spn][gene][concat_id]["seq"])
                                    otu_ = otu.replace(" ", "_")
                                    otu_ = otu_.replace(".", "").replace("'", "")
                                    if gene in self.comb_seq.keys():
                                        if otu_ not in self.comb_seq[gene].keys():
                                            self.comb_seq[gene][otu_] = seq
                                            if gene in self.comb_acc:
                                                self.comb_acc[gene][otu_] = concat_id
                                            else:
                                                self.comb_acc[gene] = {otu_: concat_id}
                                        else:
                                            self.comb_seq[gene][otu_] = seq
                                            self.comb_acc[gene][otu_] = concat_id
                                    else:
                                        self.comb_seq[gene] = {otu_: seq}
                                        self.comb_acc[gene] = {otu_: concat_id}
                                    if spn != otu:
                                        self.sp_acc_comb[spn][gene][concat_id]["new tipname"] = otu_
                                    self.otu_to_spn(spn, gene, concat_id)
                                    break
                        if spn is not None:
                            break
                if len(gene_l) == len(concat_l):
                    missing_gene = [item for item in self.genes_present if item not in gene_l]
                    for genes in missing_gene:
                        self.make_empty_seq(global_spn, genes)

    def dump(self, filename="concat_checkpoint.p"):
        """ Save a concat run as pickle.
        """
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))

    def write_otu_info(self):
        """Writes output tables to file: Makes reading important information less code heavy.

        file with all relevant GenBank info to file (otu_dict).

        :return: writes output to file
        """   
        debug("write_otu_info")
        otu_dict_keys = [
            "unique_id", "spn", "original_PS_id", "concat:status"]
        with open("{}/otu_seq_info.csv".format(self.workdir), "w") as output:
            writer = csv.writer(output)
            writer.writerow(otu_dict_keys)
            # print(self.sp_acc_comb.keys())
            for otu in self.sp_acc_comb.keys():
                # print(self.sp_acc_comb[otu].keys())
                # rowinfo = spn_writer
                for loci in self.sp_acc_comb[otu].keys():
                    # print("newinfo")
                    rowinfo = [loci]
                    # print(self.sp_acc_comb[otu][loci])
                    for concat_id in self.sp_acc_comb[otu][loci].keys():
                        rowinfo_details = deepcopy(rowinfo)
                        # print(rowinfo_details)
                        rowinfo_details.append(concat_id)
                        for item in otu_dict_keys:
                            if item in self.sp_acc_comb[otu][loci][concat_id].keys():
                                tofile = str(self.sp_acc_comb[otu][loci][concat_id][item]).replace("_", " ")
                                rowinfo_details.append(tofile)
                            else:
                                rowinfo_details.append("-")
                        # print(rowinfo_details)
                        writer.writerow(rowinfo_details)                

