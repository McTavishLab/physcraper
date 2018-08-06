#!/usr/bin/env python
import os
import subprocess
import csv
import pickle
import random

from copy import deepcopy

from dendropy import Tree, \
    DnaCharacterMatrix
from Bio import Entrez 

from __init__ import debug
from __init__ import IdDicts



class Concat(object):
    """combine several physcraper runs of the same lineage with
     different genes into one final concatenated aln and tre
     """

    def __init__(self, workdir_comb, email):
        # super(PhyscraperScrape, self).__init__()
        self.workdir = workdir_comb
        # self.gene_comb = {}
        self.sp_gi_comb = {}
        self.single_runs = {}
        self.sp_counter = {}
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.comb_seq = {}
        self.comb_gi = {}
        self.aln_all = {}
        self.aln_all_len = {}
        self.num_of_genes = 0
        self.genes_present = []
        self.tre_as_start = None
        self.short_concat_seq = None
        self.concat_tips = {}
        self.email = email

    def load_single_genes(self, workdir, pickle_fn, genename):
        """load PhyScraper class objects and make a single dict per run.
        Removes abandoned nodes first.
        """
        debug("load_single_genes: {}".format(genename))
        # debug("{}/{}".format(workdir, pickle_fn))
        scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn), 'rb'))
        scrape = self.remove_aln_tre_leaf(scrape)
        self.single_runs[genename] = deepcopy(scrape)
        return

    def remove_aln_tre_leaf(self, scrape):
        """attempt to remove all occurrences in aln and tre of otu,
         that were removed sometime in the single runs.
        """
        aln_ids = set()
        for tax in scrape.data.aln:
            aln_ids.add(tax.label)
        assert aln_ids.issubset(scrape.data.otu_dict.keys())

        treed_taxa = set()
        for leaf in scrape.data.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        # debug(treed_taxa)

        for leaf in scrape.data.tre.leaf_nodes():
            if leaf.taxon not in aln_ids:
                # debug("leaf.taxon not present in aln_ids")
                # debug(leaf.taxon)
                scrape.data.tre.prune_taxa([leaf])
                scrape.data.tre.prune_taxa_with_labels([leaf.taxon])
                scrape.data.tre.prune_taxa_with_labels([leaf])
                # scrape.data.tre.taxon_namespace.remove_taxon_label(leaf.taxon.label)

                treed_taxa.remove(leaf.taxon)
        assert treed_taxa.issubset(aln_ids)
        return scrape

    def get_taxon_info(self, key, data):
        """
        If there are no taxon information (for what ever reason) try again to obtain sp names.
        """
        # debug("get_rank_info")
        if key in data:
            if data[key] == None:

                if '^ncbi:gi' in data:
                    gi_id = data['^ncbi:gi']
                Entrez.email = self.email
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
                tax_name = data[key]
        tax_name = tax_name.replace("_", " ")
        tax_name = tax_name.replace(".", "").replace("'", "")
        tax_name = tax_name.encode('ascii')
        return tax_name

    def make_concat_id_dict(self, otu, genename, concat_id):
        """make a concat_id entry with all information
        """
        # has test
        data = self.single_runs[genename].data.otu_dict[otu]
        seq = str(self.single_runs[genename].data.aln[otu])
        if '^ot:ottTaxonName' in data:
            spn = self.get_taxon_info('^ot:ottTaxonName', data)
            if spn not in self.sp_gi_comb:
                self.sp_gi_comb[spn] = {}
            if genename not in self.sp_gi_comb[spn]:
                self.sp_gi_comb[spn][genename] = {}
            if concat_id not in self.sp_gi_comb[spn][genename]:
                # debug("make concat_id")
                if '^ncbi:gi' in data:
                    gi_id = data['^ncbi:gi']
                else:
                    gi_id = self.single_runs[genename].data.otu_dict[otu]['^ot:originalLabel']
                    # print(gi_id)
                concat_dict = {"gi_id": gi_id, "seq": seq, "spn": spn, "original_PS_id": otu,
                               "concat:status": "single run"}
                self.sp_gi_comb[spn][genename][concat_id] = concat_dict
        elif '^user:TaxonName' in data:
            spn = self.get_taxon_info('^user:TaxonName', data)
            if spn not in self.sp_gi_comb:
                self.sp_gi_comb[spn] = {}
            if genename not in self.sp_gi_comb[spn]:
                self.sp_gi_comb[spn][genename] = {}
            if concat_id not in self.sp_gi_comb[spn][genename]:
                # debug("make concat_id")
                if '^ncbi:gi' in data:
                    gi_id = data['^ncbi:gi']
                else:
                    gi_id = self.single_runs[genename].data.otu_dict[otu]['^ot:originalLabel']
                    # print(gi_id)
                concat_dict = {"gi_id": gi_id, "seq": seq, "spn": spn, "original_PS_id": otu,
                               "concat:status": "single run"}
                # debug(concat_dict)
                self.sp_gi_comb[spn][genename][concat_id] = concat_dict
        else:
            print("THERE IS A SERIOUS PROBLEM....")
        # debug(concat_dict['spn'])    
        if concat_dict['spn'] == None:
            print("THERE IS A SERIOUS PROBLEM....Number2")
            #we should never get here....
            print("spn is none")
            spn = self.get_taxon_info(gi_id)
            self.sp_gi_comb[gi_id] 
            self.sp_gi_comb[spn] =  self.sp_gi_comb[gi_id] 
            del  self.sp_gi_comb[gi_id] 
            print(some)

    def combine(self):
        """combine several PhyScraper objects to make a concatenated run dict
        """
        debug("combine")
        self.num_of_genes = len(self.single_runs)
        concat_id_counter = 1
        for genename in self.single_runs:
            self.genes_present.append(genename)
            # debug(genename)
            for otu in self.single_runs[genename].data.aln.taxon_namespace:
                concat_id = "concat_{}".format(concat_id_counter)
                self.make_concat_id_dict(otu.label, genename, concat_id)
                concat_id_counter += 1
        # debug(self.sp_gi_comb)
        return

    def sp_seq_counter(self):
        """counts how many seq per sp and genes there are
        """
        # has test
        debug("sp_seq_counter")
        # debug(self.sp_gi_comb)
        for spn in self.sp_gi_comb:
            # debug(spn)
            tmp_gene = deepcopy(self.genes_present)
            for gene in self.sp_gi_comb[spn]:
                tmp_gene.remove(gene)
                spn_new = spn.replace(" ", "_")
                if spn_new in self.sp_counter:
                    self.sp_counter[spn_new][gene] = len(self.sp_gi_comb[spn][gene])
                else:
                    self.sp_counter[spn_new] = {gene: len(self.sp_gi_comb[spn][gene])}
            for item in tmp_gene:
                if spn_new in self.sp_counter:
                    self.sp_counter[spn_new][item] = 0
                else:
                    self.sp_counter[spn_new] = {item: 0}
        debug(self.sp_counter)

    def sp_to_keep(self):
        """uses the sp_counter to make a list of sp that should be kept in concatenated alignment,
        because they are the only representative of the sp.
        """
        # has test

        debug("sp to keep")
        sp_to_keep = {}
        # debug(self.sp_counter)
        for spn in self.sp_counter:
            seq_counter = True
            not_present = 0
            for gene in self.sp_counter[spn]:
                if self.sp_counter[spn][gene] == 0:
                    # debug("sp has not all genes present")
                    seq_counter = False
                    not_present += 1
            if not seq_counter:
                sp_to_keep[spn] = not_present
        debug(sp_to_keep)
        return sp_to_keep

    def make_sp_gene_dict(self, sp_to_keep):
        """is the build around to make the dicts, that are used to make it into a dendropy aln
        """
        debug("make_sp_gene_dict")
        count = 2
        self.tmp_dict = deepcopy(self.sp_gi_comb)
        del_list = []
        while len(self.tmp_dict.keys()) >= 1:
            # debug(len(self.tmp_dict.keys()))
            # debug("start of while loop")
            del_gi = {}
            for spn in self.tmp_dict.keys():
                sp_to_keep_list = sp_to_keep.keys()
                if spn.replace(" ", "_") in sp_to_keep_list:
                    # debug("in sp_to_keep_list")
                    tmp_gene = deepcopy(self.genes_present)
                    for gene in self.tmp_dict[spn]:
                        tmp_gene.remove(gene)
                        del_gi = self.select_rnd_seq(spn, gene, del_gi, count)
                    for item in tmp_gene:
                        self.make_empty_seq(spn, item)
                    self.rm_rnd_sp(del_gi)
                    del self.tmp_dict[spn]
                else:
                    for gene in self.tmp_dict[spn]:
                        del_gi = self.select_rnd_seq(spn, gene, del_gi, count)
                    self.rm_rnd_sp(del_gi)
                self.rm_empty_spn_entries(del_gi)
        self.rename_drop_tips()


    def select_rnd_seq(self, spn, gene, del_gi, count):
        """select random seq from spn and gene to combine it with a random other one from another gene,
        but same spn
        """
        # has test

        # debug("select_rnd_seq")
        print(spn, gene, del_gi, count)
        random_gen = random.choice(list(self.tmp_dict[spn][gene]))
        # debug("random_gen")
        # debug(random_gen)
        self.sp_gi_comb[spn][gene][random_gen]["concat:status"] = "used in concat"

        seq = str(self.tmp_dict[spn][gene][random_gen]["seq"])
        # debug(self.comb_seq.keys())
        spn_ = spn.replace(" ", "_")
        spn_ = spn_.replace(".", "").replace("'", "")
        if gene in self.comb_seq.keys():
            # debug("gene in comb_seq")
            # debug(self.comb_seq[gene].keys())
            if spn_ not in self.comb_seq[gene].keys():
                # debug("sp new")
                self.comb_seq[gene][spn_] = seq
                if gene in self.comb_gi:
                    self.comb_gi[gene][spn_] = random_gen
                else:
                    self.comb_gi[gene] = {spn_: random_gen}
                if gene in del_gi.keys():
                    if spn_ not in del_gi[gene].keys():
                        del_gi[gene][spn] = random_gen
                else:
                    del_gi[gene] = {spn: random_gen}
            else:
                # debug("spn already present")
                spn_new = "{}_{}".format(spn_, count)
                while spn_new in self.comb_seq[gene].keys():

                    count += 1
                    spn_new = "{}_{}".format(spn_, count)
                    print(spn_new, count)

                self.comb_seq[gene][spn_new] = seq
                self.comb_gi[gene][spn_new] = random_gen
                self.sp_gi_comb[spn][gene][random_gen]["new tipname"] = spn_new
                if gene in del_gi.keys():
                    if spn_ not in del_gi[gene].keys():
                        del_gi[gene][spn] = random_gen
                    else:
                        del_gi[gene][spn] = random_gen
                else:
                    del_gi[gene] = {spn: random_gen}
        else:
            # debug("new gene in comb_seq")
            self.comb_seq[gene] = {spn_: seq}
            self.comb_gi[gene] = {spn_: random_gen}
            if gene in del_gi.keys():
                if spn_ not in del_gi[gene].keys():
                    del_gi[gene][spn] = random_gen
                else:
                    del_gi[gene][spn_new] = random_gen
            else:
                del_gi[gene] = {spn: random_gen}
        self.otu_to_spn(spn, gene, del_gi[gene][spn])
        debug(del_gi)
        return del_gi

    def otu_to_spn(self, spn_, gene, random_gen):
        """ makes a dict that contains the original tip labels from the starting tree, and the name it needs to have.
        This dict will be used in rename_drop_tips to rename or remove the tips. As such it is a helper function to produce the correct starting tree,
        for the concatentation run

        :param spn_: species name for concat
        :param gene: gene
        :param random_gen: the corresponding otu
        :return: self.concat_tips
        """
        # debug("otu_to_spn")
        if self.tre_start_gene == gene:
            # debug(self.tre_as_start.taxon_namespace)
            spn = spn_.replace("_", " ")
            former_otu = self.sp_gi_comb[spn][gene][random_gen]['original_PS_id']
            for otu in self.tre_as_start.taxon_namespace:
                if otu.label == former_otu:
                    if 'new tipname' in self.sp_gi_comb[spn][gene][random_gen]:
                        spn_ = self.sp_gi_comb[spn][gene][random_gen]['new tipname']
                    self.concat_tips[otu.label] = spn_
            # debug("self.concat_tips")
            # debug(self.concat_tips)
        return self.concat_tips

    def rename_drop_tips(self):
        """ remove tips from tre as start, that are not present in the concatenated aln
        rename tips that are present
        """
        debug("rename_drop_tips")
        # leaf.taxon is never in concat_tips
        for leaf in self.tre_as_start.leaf_nodes():
            if leaf.taxon.label not in self.concat_tips.keys():
                self.tre_as_start.prune_taxa([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.label])
                self.tre_as_start.prune_taxa_with_labels([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
            else:
                # debug("change label")
                for otu in self.concat_tips.keys():
                    if otu == leaf.taxon.label:
                        leaf.taxon.label = self.concat_tips[otu]

    def add_to_del_gi(self, del_gi, gene, spn, random_gen):
        """add gi number to del_gi,
        del_gi is used to remove gi's from tmp_dict, so that they will
        not be added to the concat dict twice.
        """
        spn_ = spn.replace(" ", "_")
        if gene in del_gi.keys():
            if spn_ not in del_gi[gene].keys():
                del_gi[gene][spn] = random_gen
            else:
                del_gi[gene][spn_new.format("_", " ")] = random_gen
        else:
            del_gi[gene] = {spn: random_gen}
        return del_gi

    def make_empty_seq(self, spn, gene):
        """when there is no seq for one of the genes,
        make an empty seq. Dendropy needs same taxon_namespace and number otu's for concatenation
        """
        # debug("make_empty_seq")
        for tax, seq in self.single_runs[gene].data.aln.items():
            len_gene_aln = len(seq)
            break
        empty_seq = "?" * len_gene_aln
        if gene in self.comb_seq:
            self.comb_seq[gene][spn.replace(" ", "_")] = empty_seq
        else:
            self.comb_seq[gene] = {spn.replace(" ", "_"): empty_seq}

    def rm_rnd_sp(self, del_gi):
        """removes the random selected seq from the tmp_dict
        """
        # debug("rm_rnd sp")
        for spn2 in self.tmp_dict:
            for gene2 in self.tmp_dict[spn2]:
                if gene2 in del_gi:
                    if spn2 in del_gi[gene2]:
                        key = del_gi[gene2][spn2]
                        if key in self.tmp_dict[spn2][gene2]:
                            del self.tmp_dict[spn2][gene2][key]

    def rm_empty_spn_entries(self, del_gi):
        """removes keys from tmp dict, if the key/sp has no value anymore
        """
        # debug("rm_empty_spn_entries")
        del_sp = None
        for spn2 in self.tmp_dict:
            for gene2 in self.tmp_dict[spn2]:
                if gene2 in del_gi:
                    if spn2 in del_gi[gene2]:
                        if len(self.tmp_dict[spn2][gene2]) == 0:
                            del_sp = spn2
        if del_sp is not None:
            for item in self.sp_gi_comb[del_sp]:
                for otu in self.sp_gi_comb[del_sp][item]:
                    if self.sp_gi_comb[del_sp][item][otu]["concat:status"] != "used in concat":
                        self.sp_gi_comb[del_sp][item][otu][
                            "concat:status"] = "deleted, because not enough seq are present"
            del self.tmp_dict[del_sp]

    def make_alns_dict(self):
        """make dendropy aln out of dicts for all genes
        """
        debug("make_alns_dict")
        firstelement = True
        count = 0
        for gene in self.comb_seq.keys():
            if count == 0:
                len1 = len(self.comb_seq[gene].keys())
                item_of_gene1 = self.comb_seq[gene].keys()
                item_of_gene2 = list()
                len2 = len1
                count = 1
            else:
                len2 = len(self.comb_seq[gene].keys())
                item_of_gene2 = self.comb_seq[gene].keys()
            # debug([item for item in item_of_gene1 if item not in item_of_gene2])
            # debug([item for item in item_of_gene2 if item not in item_of_gene1])
            assert len1 == len2
        for gene in self.comb_seq.keys():
            if firstelement:
                aln1 = DnaCharacterMatrix.from_dict(self.comb_seq[gene])
                firstelement = False
                self.aln_all[count] = aln1
                aln1.write(path="{}/aln_0.fas".format(self.workdir),
                           schema="fasta")
                # debug(aln1.as_string(schema="fasta"))
            else:
                aln = DnaCharacterMatrix.from_dict(self.comb_seq[gene], taxon_namespace=aln1.taxon_namespace)
                self.aln_all[count] = aln
                aln.write(path="{}/aln_{}.fas".format(self.workdir, count),
                          schema="fasta")
                # aln_all[i+1] = concat_aln
                # debug(aln.as_string(schema="fasta"))
            count += 1

    def make_concat_table(self):
        """ Make a table that shows which gi numbers were combined
        """
        genel = []
        spn_l = {}
        # debug(self.comb_gi)
        # debug(self.sp_gi_comb)
        for gene in self.comb_gi:
            genel.append(gene)
            for spn in self.comb_gi[gene]:
                concat_id = self.comb_gi[gene][spn]
                multiname = False
                if spn.split("_")[-1].isdigit():
                    tmp_spn = spn[:-2]
                    tmp_spn = tmp_spn.replace("_", " ")
                    multiname = True
                spn = spn.replace("_", " ")
                if multiname:
                    if tmp_spn in self.sp_gi_comb.keys():
                        # if self.sp_gi_comb[tmp_spn][gene][concat_id]['new tipname'] == spn:
                        gi_id = self.sp_gi_comb[tmp_spn][gene][concat_id]["gi_id"]
                    else:
                        gi_id = self.sp_gi_comb[spn][gene][concat_id]["gi_id"]
                else:
                    gi_id = self.sp_gi_comb[spn][gene][concat_id]["gi_id"]
                if spn in spn_l.keys():
                    spn_l[spn].append(gi_id)
                else:
                    spn_l[spn] = [gi_id]
        print(spn_l)
        with open('{}/concatenation.csv'.format(self.workdir), 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(genel)
            for key, value in spn_l.items():
                print(key, value)
                val_str = ''.join(str(item) for item in value)
                writer.writerow([key, val_str])

    def get_short_seq_from_concat(self, percentage = 0.37):
        """find short sequences, all below a certain threshold will be removed,
        to avoid having reallz low coverage in the aln. Default = 0.37. Note percentage is a bit misleading. 
        the cutoff is 37% of the whole concatenated alignment, but the seqeuences length is calculated without gaps present. 
        The default is so low, as I want to keep taxa that have only a single locus 
        and which is not the longest among the loci within the aln.
        """
        debug("get_short_seq_from_concat")
        seq_len = {}
        num_tax = 0
        for tax, seq in self.concatenated_aln.items():
            seq = seq.symbols_as_string().replace("-", "").replace("?", "")
            seq_len[tax] = len(seq)
            num_tax += 1

        for tax, seq in self.concatenated_aln.items():
            total_len = len(seq)
            break
        min_len = (total_len * percentage)
        prune_shortest = []
        for tax, len_seq in seq_len.items():
            if len_seq < min_len:
                prune_shortest.append(tax)
        self.short_concat_seq = prune_shortest
        # debug(self.short_concat_seq)

    def remove_short_seq(self):
        """remove short seq that were found with get_short_seq
        and write it to file
        """
        debug("remove_short_seq")
        # debug(self.tre_as_start.taxon_namespace)
        # debug(self.short_concat_seq)
        # debug(len(self.concatenated_aln))
        self.concatenated_aln.remove_sequences(self.short_concat_seq)
        # debug(len(self.concatenated_aln))
        # debug(self.tre_as_start.leaf_nodes())
        for leaf in self.tre_as_start.leaf_nodes():
            # debug(leaf.taxon.label)
            # debug(self.short_concat_seq)
            for tax in self.short_concat_seq:
                # debug(tax)
                # debug(tax.label)
                if tax.label == leaf.taxon.label.replace(" ", "_"):
                    # debug("do something")
                    self.tre_as_start.prune_taxa([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.label])
                    self.tre_as_start.prune_taxa_with_labels([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                    self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)

                else:
                    leaf.taxon.label = leaf.taxon.label.replace(" ", "_")

        debug("writing files")
        tre_as_start_str = self.tre_as_start.as_string(schema="newick",
                                                       # preserve_underscores=True,
                                                       unquoted_underscores=True,
                                                       suppress_rooting=True)

        fi = open("{}/{}".format(self.workdir, "starting_red.tre"), "w")
        fi.write(tre_as_start_str)
        fi.close()
        for tax in self.concatenated_aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")
        # debug(self.concatenated_aln.taxon_namespace)
        self.concatenated_aln.write(path="{}/{}".format(self.workdir, "concat_red.fasta"), schema="fasta")

        tre_ids = set()
        for tax in self.tre_as_start.taxon_namespace:
            tre_ids.add(tax.label)

        aln_ids = set()
        for tax in self.concatenated_aln.taxon_namespace:
            aln_ids.add(tax.label)

        # debug(len(self.otu_dict.keys()))
        # debug(len(aln_ids))
        # debug([item for item in tre_ids if item not in aln_ids])
        # debug([item for item in aln_ids if item not in tre_ids])

    def get_largest_tre(self):
        """find the single gene tree with the most tips, which will be used as
        starting tree for concat phylo reconstruction
        """
        debug('get_largest_tre')
        first = True
        len_all_taxa = {}
        for gene in self.single_runs:
            len_aln_taxa = len(self.single_runs[gene].data.aln.taxon_namespace)
            len_all_taxa[gene] = len_aln_taxa

        for gene, len_item in len_all_taxa.items():
            if first:
                len_max = len_item
                gene_max = gene
                first = False
            if len_item > len_max:
                len_max = len_item
                gene_max = gene
        self.tre_as_start = self.single_runs[gene_max].data.tre
        self.tre_start_gene = gene_max
        # debug("self.tre_as_start")
        # debug(self.tre_as_start)

    def concatenate_alns(self):
        """concatenate all alns into one aln
        """
        count = 0
        for gene in self.aln_all:
            if count == 0:
                aln1 = self.aln_all[gene]
                aln1.write(path="{}/aln1.fas".format(self.workdir),
                   schema="fasta")
                count = 1
            else:
                aln2 = self.aln_all[gene]
                count +=1
                aln2.write(path="{}/aln{}.fas".format(self.workdir, count),
                   schema="fasta")
                assert aln1.taxon_namespace == aln2.taxon_namespace
                aln1 = DnaCharacterMatrix.concatenate([aln1, aln2])
        aln1.write(path="{}/concat.fas".format(self.workdir),
                   schema="fasta")
        self.concatenated_aln = aln1

    def dump(self, filename="concat_checkpoint.p"):
        """ save a concat run as pickle
        """
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))

    def write_partition(self):
        """write the partitioning file for RAxML
        """
        debug("write_partition")
        count = 0
        for gene in self.single_runs:
            # debug(gene)
            for tax, seq in self.single_runs[gene].data.aln.items():
                len_gene = len(seq.symbols_as_string())
                break
            if count == 0:
                # debug("new file created")
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
        """places the new seqs (that are only found in loci which is not the starting tree) onto one of the single run trees
        """
        debug("place_new_seqs")
        if len(self.concatenated_aln.taxon_namespace)-len(self.short_concat_seq) > len(self.tre_as_start.leaf_nodes()):
            if os.path.exists("RAxML_labelledTree.PLACE"):
                os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
            cwd = os.getcwd()
            os.chdir(self.workdir)

            debug("make place-tree")
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-f", "v", "-q", "partition",
                             "-s", "concat_red.fasta",
                             "-t", "starting_red.tre",
                             "-n", "PLACE"])
            os.chdir(cwd)
            debug("read place tree")
            placetre = Tree.get(path="{}/starting_red.tre".format(self.workdir),
                                schema="newick",
                                preserve_underscores=True,
                                suppress_internal_node_taxa=True, suppress_leaf_node_taxa=True)
            debug("resolve polytomies")
            placetre.resolve_polytomies()
            placetre.write(path="{}/place_resolve.tre".format(self.workdir), schema="newick", unquoted_underscores=True)

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree
        """
        debug("run full tree")
        cwd = os.getcwd()
       
        os.chdir(self.workdir)
        print("after change dir")
        print(os.path.exists("concat_red.fasta.reduced"))

        if os.path.exists("place_resolve.tre"):
            starting_fn = 'place_resolve.tre'
        else:
            starting_fn = "starting_red.tre"
        if os.path.exists("concat_red.fasta.reduced"):
            aln = "concat_red.fasta.reduced"
            partition = "partition.reduced"
        else:
            aln = "concat_red.fasta"
            partition = "partition"
        
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                     "-s", aln, "--print-identical-sequences",
                     "-t", "{}".format(starting_fn),
                     "-p", "1", "-q", partition,
                     "-n", "concat"])

        if os.path.exists("concat_red.fasta.reduced"):
            aln = "concat_red.fasta.reduced"
            partition = "partition.reduced"
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                     "-s", aln, "--print-identical-sequences",
                     "-t", "{}".format(starting_fn),
                     "-p", "1", "-q", partition,
                     "-n", "concat"])

        os.chdir(cwd)

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

        print(os.getcwd())
        print(os.path.exists("place_resolve.tre"))
        if os.path.exists("place_resolve.tre"):
            starting_fn = 'place_resolve.tre'
        else:
            starting_fn = "starting_red.tre"
        print(os.path.exists("concat_red.fasta.reduced"))
        if os.path.exists("concat_red.fasta.reduced"):
            aln = "concat_red.fasta.reduced"
        else:
            aln = "concat_red.fasta"
        print(starting_fn, aln)
        debug("1")
        # run bootstrap
     
        # make bipartition tree
        # is the -f b command
        # -z specifies file with multiple trees
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", aln,
                         # "-t", "place_resolve.tre",
                         "-p", "1", "-b", "1", "-#", "autoMRE",
                         "-n", "autoMRE"])
        debug("2b")
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", aln,
                         "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                         "-n", "autoMRE"])

        # strict consensus:
        debug("4")
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-J", "STRICT",
                         "-z", "RAxML_bootstrap.autoMRE",
                         "-n", "StrictCon"])
        # majority rule:
        debug("5")
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-J", "MR",
                         "-z", "RAxML_bootstrap.all",
                         "-n", "MR"])
      