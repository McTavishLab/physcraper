#!/usr/bin/env python
import os
import subprocess

import pickle
import random

from copy import deepcopy

from dendropy import Tree, \
    DnaCharacterMatrix

from __init__ import debug


class Concat(object):
    """combine several physcraper runs of the same lineage with
     different genes into one final concatenated aln and tre"""

    def __init__(self, workdir_comb):
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
        # self.tre_as_start = None
        self.short_concat_seq = None
        self.concat_tips = {}

    def load_single_genes(self, workdir, pickle_fn, genename):
        """load Physcraper class objects and make a single dict
        """
        debug("load_single_genes: {}".format(genename))
        # debug("{}/{}".format(workdir, pickle_fn))
        scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn), 'rb'))
        scrape = self.remove_aln_tre_leaf(scrape)
        self.single_runs[genename] = deepcopy(scrape)
        return

    def remove_aln_tre_leaf(self, scrape):
        """attempt to remove all occurrences in aln and tre of otu!!
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

    def make_concat_id_dict(self, otu, genename, concat_id):
        """make a concat_id entry with all information
        """
        data = self.single_runs[genename].data.otu_dict[otu]
        seq = str(self.single_runs[genename].data.aln[otu])
        # debug(data['^ot:ottTaxonName'])
        if data['^ot:ottTaxonName'] not in self.sp_gi_comb:
            self.sp_gi_comb[data['^ot:ottTaxonName']] = {}
        if genename not in self.sp_gi_comb[data['^ot:ottTaxonName']]:
            self.sp_gi_comb[data['^ot:ottTaxonName']][genename] = {}
        if concat_id not in self.sp_gi_comb[data['^ot:ottTaxonName']][genename]:
            # debug("make concat_id")
            if '^ncbi:gi' in data:
                gi_id = data['^ncbi:gi']
            else:
                gi_id = data['^ot:ottTaxonName']
            concat_dict = {"gi_id": gi_id, "seq": seq, "spn": str(data['^ot:ottTaxonName']), "original_PS_id": otu,
                           "concat:status": "single run"}
            # debug(concat_dict)
        self.sp_gi_comb[data['^ot:ottTaxonName']][genename][concat_id] = concat_dict

    def combine(self, genelist):
        """combine several Physcraper objects to make a concatenated run dict
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
        debug(self.sp_gi_comb)
        return

    def sp_seq_counter(self):
        """counts how many seq per sp and gene there are
        """
        debug("sp_seq_counter")
        for spn in self.sp_gi_comb:
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
        because thez are the only representative of the sp
        """
        debug("sp to keep")
        # keep_sp = False
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
        # debug(sp_to_keep)
        return sp_to_keep

    def make_sp_gene_dict(self, sp_to_keep):
        """is the build around to make the dicts, that are used to make it into a dendropy aln
        """
        debug("make_sp_gene_dict")
        count = 2
        self.tmp_dict = deepcopy(self.sp_gi_comb)

        while len(self.tmp_dict.keys()) >= 1:
            print(len(self.tmp_dict.keys()))
            print("start of while loop")
            del_gi = {}

            for spn in self.tmp_dict.keys():
                debug(spn)
                # debug(sp_to_keep)
                # debug(some)
                sp_to_keep_list = sp_to_keep.keys()
                # debug(spn)
                if spn.replace(" ", "_") in sp_to_keep_list:
                    print("in sp_to_keep_list")
                    tmp_gene = deepcopy(self.genes_present)
                    for gene in self.tmp_dict[spn]:
                        debug(gene)
                        tmp_gene.remove(gene)
                        # if gene in self.tmp_dict[spn]:
                        #     if len(self.tmp_dict[spn][gene]) >= 1:
                        print("if gene has seq...")
                        # debug(gene)
                        del_gi = self.select_rnd_seq(spn, gene, del_gi, count)
                    for item in tmp_gene:
                        print('tmp_gene')

                        print(tmp_gene)
                        # print("make empty seq")
                        print(item)

                        self.make_empty_seq(spn, item)
                    debug(del_gi)
                    self.rm_rnd_sp(del_gi)
                    del self.tmp_dict[spn]

                else:
                    print("im in else")
                    for gene in self.tmp_dict[spn]:
                        del_gi = self.select_rnd_seq(spn, gene, del_gi, count)

                # print("comb)")
                # print(self.comb_gi)
                # # print(self.comb_seq)
                # # # del_gi = comb_gi
                # print("del gi")
                # print(del_gi)

                self.rm_rnd_sp(del_gi)
                # #remove empty spn entries frmo tmp_d
                # for spn2 in self.tmp_dict:
                #     print(spn2)
                #     for gene2 in self.tmp_dict[spn2]:
                #         print(gene2)
                #         # print()
                #         # for key in comb_gi[gene][spn]:
                #         #     # print(key)
                #         #     # print(self.tmp_dict[spn][gene])
                #         spn2_ = spn2.replace(" ","_")
                #         print(del_gi[gene2])
                #         if spn2 in del_gi[gene2]:
                #             key = del_gi[gene2][spn2]
                #             print(key)
                #             if key in self.tmp_dict[spn2][gene2]:
                #                 print("del")
                #                 del self.tmp_dict[spn2][gene2][key]
                #                 if len(self.tmp_dict[spn2][gene2]) == 0:
                #                     # print("nested dict is empty")
                #                     del_list.append(spn2)

                # # print("self.tmp_dict")
                # # print(self.tmp_dict.keys())
                # # print("delllist")
                # del_list = set(del_list)
                # # print(del_list)
                # for sp_to_del in del_list:
                #     # print(sp_to_del)
                #     del self.tmp_dict[sp_to_del]
                # print(self.tmp_dict)
                # print(len(self.tmp_dict.keys()))
                self.otu_to_spn(spn, gene, del_gi[gene][spn])
                self.rm_empty_spn_entries(del_gi)
                print("end of while loop")
        debug(self.comb_seq)

        self.remove_drop_tips()

        # print(some)

    # def make_sp_gene_dict(self, sp_to_keep):
    #     """is the build around to make the dicts,
    #     that are used to make it into a concatenated dendropy aln
    #     """
    #     debug("make_sp_gene_dict")
    #     count = 2
    #     self.tmp_dict = self.sp_gi_comb
    #     print(self.sp_gi_comb)
    #     self.concat_tips = {}
    #     while len(self.tmp_dict.keys()) >= 1:
    #         # print(len(self.tmp_dict.keys()))
    #         # print("start of while loop")
    #         del_gi = {}

    #         for spn in self.tmp_dict.keys():
    #             debug(spn)
    #             # debug(sp_to_keep)
    #             # debug(some)
    #             sp_to_keep_list = sp_to_keep.keys()
    #             debug("sp_to_keep_list")
    #             debug(sp_to_keep_list)

    #             if spn.replace(" ", "_") in sp_to_keep_list:
    #                 # print("do something")
    #                 tmp_gene = deepcopy(self.genes_present)
    #                 for gene in self.tmp_dict[spn]:
    #                     tmp_gene.remove(gene)
    #                     # if gene in self.tmp_dict[spn]:
    #                 #     if len(self.tmp_dict[spn][gene]) >= 1:
    #                     # print("if gene has seq...")
    #                     # debug(gene)
    #                     del_gi = self.select_rnd_seq(spn, gene, del_gi, count)
    #                 for item in tmp_gene:
    #                     # print("make empty seq")
    #                     # print(item)

    #                     self.make_empty_seq(spn, item)
    #                 # debug(del_gi)
    #                 self.rm_empty_spn_entries(del_gi)
    #             else:
    #                 print("im in else make_sp_gene_dict")
    #                 for gene in self.tmp_dict[spn]:
    #                     del_gi = self.select_rnd_seq(spn, gene, del_gi, count)

    #             self.rm_empty_spn_entries(del_gi)
    #             print("end of while loop")
    #     debug(del_gi)
    #     # self.remove_drop_tips()
    #     # print("self.tre_as_start.taxon_namespace")
    #     # print(self.tre_as_start.taxon_namespace)
    #     # debug(some)
    #     print(self.comb_seq)
    #     # debug(some)

    def select_rnd_seq(self, spn, gene, del_gi, count):
        """select random seq from spn and gene to combine it with a random other one from another gene,
        but same spn
        """
        debug("select_rnd_seq")
        random_gen = random.choice(list(self.tmp_dict[spn][gene]))
        print("random_gen")
        print(random_gen)
        self.sp_gi_comb[spn][gene][random_gen]["concat:status"] = "used in concat"

        seq = str(self.tmp_dict[spn][gene][random_gen]["seq"])
        print(self.comb_seq.keys())
        spn_ = spn.replace(" ", "_")
        if gene in self.comb_seq.keys():
            print("gene in comb_seq")
            print(spn_, gene)
            print(self.comb_seq[gene].keys())
            if spn_ not in self.comb_seq[gene].keys():
                print("sp new")
                self.comb_seq[gene][spn_] = seq
                if gene in self.comb_gi:
                    self.comb_gi[gene][spn_] = random_gen
                else:
                    self.comb_gi[gene] = {spn_: random_gen}
                if gene in del_gi.keys():
                    if spn_ not in del_gi[gene].keys():
                        del_gi[gene][spn] = random_gen
                    else:
                        del_gi[gene][spn_new.format("_", " ")] = random_gen
                else:
                    del_gi[gene] = {spn: random_gen}
            else:
                print("spn presemt")
                spn_new = "{}_{}".format(spn_, count)
                while spn_new in self.comb_seq[gene].keys():
                    print("spn already present")
                    count += 1
                    spn_new = "{}_{}".format(spn_, count)

                self.comb_seq[gene][spn_new] = seq
                self.comb_gi[gene][spn_new] = random_gen
                if gene in del_gi.keys():
                    if spn_ not in del_gi[gene].keys():
                        del_gi[gene][spn] = random_gen
                    else:
                        del_gi[gene][spn_new] = random_gen
                else:
                    del_gi[gene] = {spn: random_gen}
        else:
            print("new gene in comb_seq")
            self.comb_seq[gene] = {spn_: seq}
            self.comb_gi[gene] = {spn_: random_gen}
            if gene in del_gi.keys():
                if spn_ not in del_gi[gene].keys():
                    del_gi[gene][spn] = random_gen
                else:
                    del_gi[gene][spn_new] = random_gen
            else:
                del_gi[gene] = {spn: random_gen}
        return del_gi

    # def select_rnd_seq(self, spn, gene, del_gi, count):
    #     """select random seq from spn and gene to combine it with a random other one from another gene, but same spn
    #     """
    #     debug("select_rnd_seq")
    #     # print(self)
    #     # print(type(self))
    #     random_gen = random.choice(list(self.tmp_dict[spn][gene]))
    #     self.sp_gi_comb[spn][gene][random_gen]["concat:status"] = "used in concat"

    #     seq = str(self.tmp_dict[spn][gene][random_gen]["seq"])
    #     # print(comb_seq.keys())
    #     spn_ = spn.replace(" ", "_")
    #     if gene in self.comb_seq.keys():
    #         print(spn_, gene)
    #         print(self.comb_seq[gene].keys())
    #         if spn_ not in self.comb_seq[gene].keys():
    #             print("in if")
    #             self.comb_seq[gene][spn_] = seq
    #             if gene in self.comb_gi:
    #                 self.comb_gi[gene][spn_] = random_gen
    #             else:
    #                 self.comb_gi[gene] ={spn_: random_gen}
    #             del_gi = self.add_to_del_gi(del_gi, gene, spn, random_gen)
    #         else:
    #             print("in else")
    #             spn_new ="{}_{}".format(spn_, count)
    #             while spn_new in self.comb_seq[gene].keys():
    #                 count += 1
    #                 spn_new ="{}_{}".format(spn_, count)
    #                 debug(spn_new)
    #             # spn_ = spn_new
    #             print(spn_new)

    #             self.comb_seq[gene][spn_new] = seq
    #             self.comb_gi[gene][spn_new] = random_gen
    #             del_gi = self.add_to_del_gi(del_gi, gene, spn_new, random_gen)
    #     else:
    #         self.comb_seq[gene] = {spn_: seq}
    #         self.comb_gi[gene] = {spn_: random_gen}
    #         del_gi = self.add_to_del_gi(del_gi, gene, spn, random_gen)
    #     # print(self.comb_seq)

    #     ### transform the tipnames of the concat starting tree

    #             # spn_not_present = False
    #         # if spn_not_present == True:
    #         #     for leaf in self.tre_as_start.leaf_nodes():
    #         #         print(leaf, leaf.taxon)
    #         #         debug(some)
    #         #         # if leaf == otu:
    #             # self.tre_as_start.prune_taxa([leaf])
    #             # self.tre_as_start.prune_taxa_with_labels([leaf.taxon])
    #             # self.tre_as_start.prune_taxa_with_labels([leaf])

    #     print(spn_, random_gen)
    #     self.concat_tips = self.otu_to_spn(spn_, gene, random_gen)
    #     self.sp_gi_comb[spn][gene][random_gen]

    #     return del_gi

    def otu_to_spn(self, spn_, gene, random_gen):
        debug("otu_to_spn")
        print(spn_, gene)
        # print(self.tre_as_start.taxon_namespace)
        spn = spn_.replace("_", " ")
        former_otu = self.sp_gi_comb[spn][gene][random_gen]['original_PS_id']
        # debug("former_otu")
        # debug(former_otu)
        # spn_not_present = True
        # print(type(self.tre_as_start))
        for otu in self.tre_as_start.taxon_namespace:
            # debug(otu)
            if otu.label == former_otu:
                # debug("names are equal")
                new_otu = spn_
                # otu.label = spn_
                self.concat_tips[otu.label] = spn_
        # print("self.concat_tips")
        print(self.concat_tips)

        return self.concat_tips

    def remove_drop_tips(self):

        debug("remove_drop_tips")
        debug(self.concat_tips)
        ## leaf.taxon is never in concat_tips
        for leaf in self.tre_as_start.leaf_nodes():
            # debug(leaf.taxon)
            # debug(self.concat_tips.keys())
            if leaf.taxon.label not in self.concat_tips.keys():
                debug("prune")

                # debug("leaf.taxon not present in aln_ids")
                # debug(leaf.taxon)
                self.tre_as_start.prune_taxa([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.label])
                self.tre_as_start.prune_taxa_with_labels([leaf])
                self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
            else:
                # debug("change label")
                for otu in self.concat_tips.keys():
                    if otu == leaf.taxon.label:
                        debug("equal")
                        leaf.taxon.label = self.concat_tips[otu]

        # debug(self.tre_as_start.taxon_namespace)
        # debug(self.tre_as_start)
        # debug(some)

        # remove_l = []
        # for otu in self.tre_as_start.taxon_namespace:
        #     # debug(otu)

        #     # debug(otu.label)
        #     if otu.label in self.concat_tips.keys():
        #         debug("change label")
        #         otu.label = self.concat_tips[otu.label]
        #     else:
        #         remove_l.append(otu.label)
        #         debug("prune")
        #         debug(otu)
        #         debug(otu.label)
        #         print(type(self.tre_as_start))

        #         self.tre_as_start.prune_taxa([otu])
        #         self.tre_as_start.prune_taxa_with_labels([otu.label])
        #         self.tre_as_start.prune_taxa_with_labels([otu])
        # self.tre_as_start.prune_taxa(remove_l)
        # self.tre_as_start.prune_taxa_with_labels(remove_l)
        # self.tre_as_start.prune_taxa_with_labels([otu.label])

        debug(self.tre_as_start.taxon_namespace)
        # debug(seom)

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
        debug("make_empty_seq")
        for tax, seq in self.single_runs[gene].data.aln.items():
            len_gene_aln = len(seq)
            break
        empty_seq = "?" * len_gene_aln
        if gene in self.comb_seq:
            self.comb_seq[gene][spn.replace(" ", "_")] = empty_seq
        else:
            self.comb_seq[gene] = {spn.replace(" ", "_"): empty_seq}

    def rm_rnd_sp(self, del_gi):
        """removes keys from tmp dict, if the key/sp has no value anymore
        """
        debug("rm_empty_spn_entries")
        debug(del_gi)
        for spn2 in self.tmp_dict:
            # print(spn2)
            for gene2 in self.tmp_dict[spn2]:
                # print(gene2)
                # print()
                # for key in comb_gi[gene][spn]:
                #     # print(key)
                #     # print(self.tmp_dict[spn][gene])
                spn2_ = spn2.replace(" ", "_")

                if gene2 in del_gi:
                    # print(del_gi[gene2])
                    if spn2 in del_gi[gene2]:
                        key = del_gi[gene2][spn2]
                        print(key)
                        if key in self.tmp_dict[spn2][gene2]:
                            print("del")
                            print(self.tmp_dict[spn2][gene2][key])
                            del self.tmp_dict[spn2][gene2][key]
                            # for genes in self.tmp_dict[spn2]:
                            #     if len(self.tmp_dict[spn2][genes]) == 0:
                            #     # print("nested dict is empty")
                            #         del_list.append(spn2)

    def rm_empty_spn_entries(self, del_gi):
        """removes keys from tmp dict, if the key/sp has no value anymore
        """
        debug("rm_empty_spn_entries")
        del_list = []
        debug(del_gi)
        # for gene in del_gi.keys():
        #     for spn in del_gi[gene]:
        #         for rnd_seq  in del_gi[gene][spn]:
        #             debug(self.tmp_dict[spn][gene])
        #             if len(self.tmp_dict[spn][gene]) == 0:
        #                 del_list.append(spn)

        for spn2 in self.tmp_dict:
            # debug(spn2)
            for gene2 in self.tmp_dict[spn2]:
                # debug(gene2)
                # for key in comb_gi[gene][spn]:
                #     # debug(key)
                #     # debug(self.tmp_dict[spn][gene])
                spn2_ = spn2.replace(" ", "_")

                if gene2 in del_gi:
                    # debug("gene2 in del_gi")
                    # debug(del_gi[gene2])
                    if spn2 in del_gi[gene2]:
                        # debug("spn2 in del_gi")
                        key = del_gi[gene2][spn2]
                        debug(self.tmp_dict[spn2][gene2])
                        if len(self.tmp_dict[spn2][gene2]) == 0:
                            del_list.append(spn2)

        # print("self.tmp_dict")
        # print(self.tmp_dict.keys())
        # print("delllist")
        del_list = set(del_list)
        # print(del_list)

        for sp_to_del in del_list:
            print(sp_to_del)

            for item in self.sp_gi_comb[sp_to_del]:
                # for gene in self.sp_gi_comb[sp_to_del]:
                for otu in self.sp_gi_comb[sp_to_del][item]:
                    if self.sp_gi_comb[sp_to_del][item][otu]["concat:status"] != "used in concat":
                        # debug("change status")
                        debug(self.sp_gi_comb[sp_to_del][item][otu])
                        self.sp_gi_comb[sp_to_del][item][otu][
                            "concat:status"] = "deleted, because not enough seq are present"

            del self.tmp_dict[sp_to_del]
            if sp_to_del == "Senecio cf. lautus DOB-2013":
                debug(self.comb_seq)
        debug(self.tmp_dict.keys())

    # def rm_empty_spn_entries(self, del_gi):
    #     """removes keys from tmp dict, if the key/sp has no value anymore
    #     """
    #     debug("rm_empty_spn_entries")
    #     del_list = []
    #     debug(del_gi)
    #     for spn2 in self.tmp_dict:
    #         for gene2 in self.tmp_dict[spn2]:
    #             spn2_ = spn2.replace(" ","_")
    #             if gene2 in del_gi:
    #                 if spn2 in del_gi[gene2]:
    #                     key = del_gi[gene2][spn2]
    #                     if key in self.tmp_dict[spn2][gene2]:
    #                         # print("del")
    #                         del self.tmp_dict[spn2][gene2][key]
    #                         if len(self.tmp_dict[spn2][gene2]) == 0:
    #                             print("nested dict is empty")
    #                         del_list.append(spn2)
    #     del_list = set(del_list)
    #     for sp_to_del in del_list:
    #         print(sp_to_del)
    #         for item in self.sp_gi_comb[sp_to_del]:
    #             for gene in self.sp_gi_comb[sp_to_del]:
    #                 for otu in self.sp_gi_comb[sp_to_del][item]:
    #                     if self.sp_gi_comb[sp_to_del][item][otu]["concat:status"] != "used in concat":
    #                         # debug("change status")
    #                         self.sp_gi_comb[sp_to_del][item][otu]["concat:status"] = "deleted, because not enough seq are present"

    #         # self.sp_gi_comb[spn2][gene2]
    #         del self.tmp_dict[sp_to_del]
    #     # debug(self.tmp_dict.keys())
    #     debug("del_list")
    #     debug(del_list)
    #     # some

    def make_alns_dict(self):
        """make dendropy aln out of dicts for all genes
        """
        debug("make_alns_dict")
        count = 1
        firstelement = False
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
            debug([item for item in item_of_gene1 if item not in item_of_gene2])
            debug([item for item in item_of_gene2 if item not in item_of_gene1])
            assert len1 == len2
        for gene in self.comb_seq.keys():
            if not firstelement:
                aln1 = DnaCharacterMatrix.from_dict(self.comb_seq[gene])
                firstelement = True
                self.aln_all[count] = aln1
                aln1.write(path="{}/aln_0.fas".format(self.workdir),
                           schema="fasta")
                # print(aln1.as_string(schema="fasta"))
            else:
                aln = DnaCharacterMatrix.from_dict(self.comb_seq[gene], taxon_namespace=aln1.taxon_namespace)
                self.aln_all[count] = aln
                aln.write(path="{}/aln_{}.fas".format(self.workdir, count),
                          schema="fasta")

                # aln_all[i+1] = concat_aln
                print(aln.as_string(schema="fasta"))
            count += 1

    def get_short_seq_from_concat(self):
        """find short sequences, all below a certain threshold will be removed,
        to avoid having reallz low coverage in the aln
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

        min_len = (total_len * 0.37)

        prune_shortest = []
        for tax, len_seq in seq_len.items():
            if len_seq < min_len:
                prune_shortest.append(tax)
        self.short_concat_seq = prune_shortest
        debug(self.short_concat_seq)

    def remove_short_seq(self):
        """remove short seq that were found with get_short_seq
        """
        debug("remove_short_seq")
        debug(self.short_concat_seq)
        debug(len(self.concatenated_aln))
        self.concatenated_aln.remove_sequences(self.short_concat_seq)
        debug(len(self.concatenated_aln))
        debug(self.tre_as_start.leaf_nodes())
        for leaf in self.tre_as_start.leaf_nodes():
            debug(leaf.taxon.label)
            debug(self.short_concat_seq)
            for tax in self.short_concat_seq:
                debug(tax)
                debug(tax.label)
                if tax.label == leaf.taxon.label.replace(" ", "_"):
                    print("do something")
                    self.tre_as_start.prune_taxa([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.label])
                    self.tre_as_start.prune_taxa_with_labels([leaf])
                    self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
                    self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
            # if leaf.taxon in self.short_concat_seq:
            #     debug("prune")
            #
            #     # debug("leaf.taxon not present in aln_ids")
            #     # debug(leaf.taxon)
            #     self.tre_as_start.prune_taxa([leaf])
            #     self.tre_as_start.prune_taxa_with_labels([leaf.label])
            #     self.tre_as_start.prune_taxa_with_labels([leaf])
            #     self.tre_as_start.prune_taxa_with_labels([leaf.taxon.label])
            #     self.tre_as_start.taxon_namespace.remove_taxon_label(leaf.taxon.label)
                else:
                    leaf.taxon.label = leaf.taxon.label.replace(" ", "_")

        print("writing files")
        tre_as_start_str = self.tre_as_start.as_string(schema="newick",
                                                       # preserve_underscores=True,
                                                       unquoted_underscores=True,
                                                       suppress_rooting=True)

        fi = open("{}/{}".format(self.workdir, "starting_red.tre"), "w")
        fi.write(tre_as_start_str)
        fi.close()
        for tax in self.concatenated_aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")
        debug(self.concatenated_aln.taxon_namespace)
        self.concatenated_aln.write(path="{}/{}".format(self.workdir, "concat_red.fasta"), schema="fasta")

        tre_ids = set()
        for tax in self.tre_as_start.taxon_namespace:
            tre_ids.add(tax.label)

        aln_ids = set()
        for tax in self.concatenated_aln.taxon_namespace:
            aln_ids.add(tax.label)

        # debug(len(self.otu_dict.keys()))
        # debug(len(aln_ids))
        debug([item for item in tre_ids if item not in aln_ids])
        debug([item for item in aln_ids if item not in tre_ids])

    def get_largest_tre(self):
        """find the single gene tree with the most tips, which will be used as
        starting tree for concat phzlo reconstruction
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
        debug("self.tre_as_start")
        debug(self.tre_as_start)

    def concatenate_alns(self):
        """concatenate all alns into one aln
        """
        # aln1_taxlist = []
        # aln2_taxlist = []
        count = 0
        for gene in self.aln_all:
            debug(gene)
            if count == 0:
                aln1 = self.aln_all[gene]
                count = 1
            else:
                aln2 = self.aln_all[gene]
                assert aln1.taxon_namespace == aln2.taxon_namespace
                # for tax in aln1.taxon_namespace:
                #     aln1_taxlist.append(tax.label)
                # for tax in aln2.taxon_namespace:
                #     aln2_taxlist.append(tax.label)
                # assert aln1_taxlist == aln2_taxlist
                aln1 = DnaCharacterMatrix.concatenate([aln1, aln2])
        aln1.write(path="{}/concat.fas".format(self.workdir),
                   schema="fasta")
        # debug(aln1.as_string(schema="fasta"))
        self.concatenated_aln = aln1

    def dump(self, filename="concat_checkpoint.p"):
        """ save a concat run as pickle
        """
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))

    def write_partition(self, gene, count):
        """write the partitioning file for RAxML
        """
        if count == 0:
            partition = open("{}/partition".format(self.workdir), "w")
            partition.close()
            with open("{}/partition".format(self.workdir), "a") as partition:
                partition.write("DNA, {} = 1-{}\n".format(gene, self.aln_all_len[gene]))
            # debug("len of 1. gnee")
            # debug(self.aln_all_len[gene])
            self.part_len = self.aln_all_len[gene]
        else:
            start = self.part_len + 1
            # debug(self.part_len, gene)
            # debug(self.aln_all_len.keys())
            end = self.part_len + self.aln_all_len[gene]
            self.part_len = self.part_len + self.aln_all_len[gene]
            with open("{}/partition".format(self.workdir), "a") as partition:
                partition.write("DNA, {} = {}-{}\n".format(gene, start, end))

    def place_new_seqs(self):
        """places the new seqs onto one of the single run trees
        """
        cwd = os.getcwd()
        debug(cwd)
        os.chdir(self.workdir)
        debug(self.workdir)
        debug("make place-tree")
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-f", "v", "-q", "partition",
                         "-s", "concat_red.fasta",
                         "-t", "starting_red.tre",
                         "-n", "PLACE"])
        # debug(some)
        os.chdir(cwd)

        debug("read place tree")
        placetre = Tree.get(path="{}/starting.tre".format(self.workdir),
                            schema="newick",
                            preserve_underscores=True)
        debug("resolve polytomies")
        placetre.resolve_polytomies()
        placetre.write(path="{}/place_resolve.tre".format(self.workdir), schema="newick", unquoted_underscores=True)
        self._query_seqs_placed = 1
        # cwd = os.getcwd()

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree"""
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                         "-s", "{}/concat_red.fasta".format(self.workdir),
                         "-t", "{}/place_resolve.tre".format(self.workdir),
                         "-p", "1", "-q", "{}/partition".format(self.workdir),
                         "-n", "concat"])
