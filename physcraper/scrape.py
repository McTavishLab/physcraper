import sys
import re
import os
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import contextlib
import time
import csv

from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel

from Bio.Blast import NCBIXML

from copy import deepcopy
import physcraper.AWSWWW as AWSWWW

from physcraper import ncbi_data_parser, filter_by_local_blast
from physcraper.configobj import ConfigObj
from physcraper.ids import IdDicts
from physcraper.aligntreetax import AlignTreeTax
from physcraper.helpers import cd


from . import writeinfofiles

_VERBOSE = 1
_DEBUG = 1
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)



class PhyscraperScrape(object):
    """
    This is the class that does the perpetual updating

        To build the class the following is needed:

          * **data_obj**: Object of class ATT (see above)
          * **ids_obj**: Object of class IdDict (see above)

        During the initializing process the following self.objects are generated:

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
          * **self.mrca_ncbi**:  int ncbi identifier of mrca

          * **self.tmpfi**: path to a file or folder???
          * **self.blast_subdir**: path to folder that contains the files writen during blast

          * **self.newseqs_file**: filename of files that contains the sequences from self.new_seqs_otu_id
          * **self.date**: Date of the run - may lag behind real date!
          * **self.repeat**: either 1 or 0, it is used to determine if we continue updating the tree, no new seqs found = 0
          * **self.newseqs_acc**: list of all gi_ids that were passed into remove_identical_seq(). Used to speed up adding process
          * **self.blacklist**: list of gi_id of sequences that shall not be added or need to be removed. Supplied by user.
          * **self.seq_filter**: list of words that may occur in otu_dict.status and which shall not be used in the building of FilterBlast.sp_d (that's the main function), but it is also used as assert statement to make sure unwanted seqs are not added.
          * **self.unpublished**: True/False. Used to look for local unpublished seq that shall be added if True.
          * **self.path_to_local_seq:** Usually False, contains path to unpublished sequences if option is used.

        Following functions are called during the init-process:

            * **self.reset_markers()**: adds things to self: I think they are used to make sure certain function run, if program crashed and pickle file is read in.
                * self._blasted: 0/1, if run_blast_wrapper() was called, it is set to 1 for the round.
                * self._blast_read: 0/1, if read_blast_wrapper() was called, it is set to 1 for the round.
                * self._identical_removed: 0
                * self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
                * self._query_seqs_aligned: 0
                * self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
                * self._reconciled: 0
                * self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
            * **self.OToL_unmapped_tips()**: function that either removes or maps unmapped taxa from OToL studies
    """
    def __init__(self, data_obj, ids_obj, ingroup_mrca=None):
        assert isinstance(data_obj, AlignTreeTax)
        assert isinstance(ids_obj, IdDicts)
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = self.ids.config  # pointer to config
        self.new_seqs = {}  # all new seq after read_blast_wrapper
        self.new_seqs_otu_id = {}  # only new seq which passed remove_identical
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)  # TODO: For what do we want to use this? Unused!
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = ""
        self.date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        self.repeat = 1  # used to determine if we continue updating the tree
        self.newseqs_acc = []  # all ever added Genbank accession numbers during any PhyScraper run, used to speed up adding process
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,",
                           "local"]  # TODO MK: try to move completely to FilterBlast class
        self.reset_markers()
        self.unpublished = False  # used to look for local unpublished seq that shall be added.
        self.path_to_local_seq = False  # path to unpublished seq.
        self.backbone = False
        self.OToL_unmapped_tips()  # added to do stuff with un-mapped tips from OToL #WTF
        self.gb_not_added = []  # list of blast seqs not added
        self.del_superseq = set()  # items that were deleted bc they are superseqs, needed for assert statement
        self.mrca_ott = data_obj.mrca_ott
        self.mrca_ncbi = self.ids.ott_to_ncbi.get(data_obj.mrca_ott)
        if self.mrca_ncbi == None:
            sys.stderr.write("ingroup mrca{} does not have a direct match to ncbi.\n".format(ingroup_mrca))
        debug("created physcraper ncbi_mrca {},".format(self.mrca_ncbi))
        self.map_taxa_to_ncbi()
        assert self.mrca_ncbi
#markers for status
#holders for new data
        self.blacklist = []
                     
    def map_taxa_to_ncbi(self):
        for otu in self.data.otu_dict:
            if self.data.otu_dict[otu].get("^ncbi:taxon") == None:
                if self.data.otu_dict[otu].get("^ot:ottId"):
                    ottid = self.data.otu_dict[otu]["^ot:ottId"]
                    self.data.otu_dict[otu]["^ncbi:taxon"]=self.ids.ott_to_ncbi.get(ottid,0)


    # TODO is this the right place for this? MK: According to PEP8, no...
    def reset_markers(self):
        self._blasted = 0
        self._blast_read = 0
        self._query_seqs_written = 0
        self._query_seqs_placed = 0
        self._full_tree_est = 0

    def reset_new_seqs_acc(self):
        """ Needs to be reseted if you want to rerun the filtering to get lower rank taxa added"""
        self.newseqs_acc = []

    def OToL_unmapped_tips(self):
        """Assign names or remove tips from aln and tre that were not mapped during initiation of ATT class.
        """
        debug("OTOL unmapped")
        if self.config.unmapped == "remove":
            for key in self.data.otu_dict:
                if "^ot:ottId" not in self.data.otu_dict[key]:
                    # second condition for OToL unmapped taxa, not present in own_data
                    if u"^ot:treebaseOTUId" in self.data.otu_dict[key]:
                        self.data.remove_taxa_aln_tre(key)
        elif self.config.unmapped == "keep":
            i = 1
            for key in self.data.otu_dict:
                i = i + 1
                if "^ot:ottId" not in self.data.otu_dict[key]:
                    self.data.otu_dict[key]["^ot:ottId"] = self.data.mrca_ott
                    if self.data.mrca_ott in self.ids.ott_to_name:
                        self.data.otu_dict[key]['^ot:ottTaxonName'] = self.ids.ott_to_name[self.data.mrca_ott]
                    else:
                        debug("think about a way...")
                        tx = APIWrapper().taxomachine
                        nms = tx.taxon(self.data.mrca_ott)
                        taxon_name = nms[u'unique_name']
                        self.data.otu_dict[key]['^ot:ottTaxonName'] = "unknown_{}".format(taxon_name)

    def run_local_blast_cmd(self, query, taxon_label, fn_path):
        """Contains the cmds used to run a local blast query, which is different from the web-queries.

        :param query: query sequence
        :param taxon_label: corresponding taxon name for query sequence
        :param fn_path: path to output file for blast query result

        :return: runs local blast query and writes it to file
        """
        abs_blastdir = os.path.abspath(self.blast_subdir)
        abs_fn = os.path.abspath(fn_path)
        toblast = open("{}/tmp.fas".format(os.path.abspath(self.blast_subdir)), "w+")
        toblast.write(">{}\n".format(taxon_label))
        toblast.write("{}\n".format(query))
        toblast.close()
        assert os.path.isdir(self.config.blastdb), ("blast dir does not exist: '%s'." % self.config.blastdb)
        with cd(self.config.blastdb):
            # this format (6) allows to get the taxonomic information at the same time
            outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid'"
            # outfmt = " -outfmt 5"  # format for xml file type
            # TODO query via stdin
            blastcmd = "blastn -query " + "{}/tmp.fas".format(abs_blastdir) + \
                       " -db {}nt -out ".format(self.config.blastdb) + abs_fn + \
                       " {} -num_threads {}".format(outfmt, self.config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(self.config.hitlist_size,
                                                                  self.config.hitlist_size)
            os.system(blastcmd)

    def local_blast_for_unpublished(self, query, taxon):
        """
        Run a local blast search if the data is unpublished.

        :param query: query sequence
        :param taxon: taxon.label used as identifier for the sequences
        :return: xml files with the results of the local blast
        """
        with cd(os.path.join(self.workdir, "blast")):
            debug("run against local unpublished data")
            debug(self.blast_subdir)
            toblast = open("{}/tmp.fas".format(self.blast_subdir), "w")
            toblast.write(">{}\n".format(taxon))
            toblast.write("{}\n".format(query))
            toblast.close()
            blast_db = "local_unpubl_seq_db"
            output = "tst_fn"
            blastcmd = "blastn -query {}/tmp.fas -db {} -out output_{}.xml " \
                       "-outfmt 5".format(self.blast_subdir, blast_db, output)
            os.system(blastcmd)

    def run_web_blast_query(self, query, equery, fn_path):
        """Equivalent to run_local_blast_cmd() but for webqueries,
        that need to be implemented differently.

        :param query: query sequence
        :param equery: method to limit blast query to mrca
        :param fn_path: path to output file for blast query result
        :return: runs web blast query and writes it to file
        """
        if self.config.url_base:
            debug("blasting {} using {}".format(self.config.url_base))
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          url_base=self.config.url_base,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size,
                                          num_threads=self.config.num_threads)
        else:
            debug("blasting {} using webservice".format(fn_path))
            result_handle = AWSWWW.qblast("blastn",
                                          "nt",
                                          query,
                                          entrez_query=equery,
                                          hitlist_size=self.config.hitlist_size)
        save_file = open(fn_path, "w")
        save_file.write(result_handle.read())
        result_handle.close()
        save_file.close()

    def run_blast_wrapper(self):  # TODO Should this be happening elsewhere?
        """generates the blast queries and saves them depending on the blasting method to different file formats

        It runs blast if the sequences was not blasted since the user defined threshold in the config file (delay).

        :return: writes blast queries to file
        """
        delay = self.config.delay
        debug("run_blast_wrapper")
        debug(self.blast_subdir)
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            if otu_id in self.data.otu_dict:
                if _VERBOSE:
                    sys.stdout.write("blasting {}\n".format(otu_id))
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(
                    last_blast, "%Y/%m/%d")).days)
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                if self.unpublished:
                    self.local_blast_for_unpublished(query, taxon.label)
                    if self.backbone is True:
                        self.data.otu_dict[otu_id]["^physcraper:last_blasted"] = today
                else:
                    if time_passed > delay:
                        if self.config.blast_loc == "local":
                            file_ending = "txt"
                        else:
                            file_ending = "xml"
                        if self.config.gb_id_filename is True:
                            fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label)
                            fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                        else:
                            fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                        # if _DEBUG:
                        #     sys.stdout.write("attempting to write {}\n".format(fn_path))
                        if not os.path.isfile(fn_path):
                            if _VERBOSE:
                                sys.stdout.write("blasting seq {}\n".format(taxon.label))
                            if self.config.blast_loc == 'local':
                                self.run_local_blast_cmd(query, taxon.label, fn_path)
                            if self.config.blast_loc == 'remote':
                                equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi, last_blast, today)
                                self.run_web_blast_query(query, equery, fn_path)
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        else:
                            if _DEBUG:
                                sys.stdout.write("file {} exists in current blast run. Will not blast, "
                                                 "delete file to force\n".format(fn_path))
                    else:
                        if _VERBOSE:
                            sys.stdout.write("otu {} was last blasted {} days ago and is not being re-blasted. "
                                             "Use run_blast_wrapper(delay = 0) to force a search.\n".format(otu_id,
                                                                                                            last_blast))
        self._blasted = 1

    def read_local_blast_query(self, fn_path):
        """ Implementation to read in results of local blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gb_dict dictionaries
        """
        # debug("read_local_blast_query")
        query_dict = {}
        with open(fn_path, mode="r") as infile:
            for lin in infile:
                # debug("new lin")
                # sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, stitle = lin.strip().split('\t')
                sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, salltitles, sallseqid = lin.strip().split('\t')
                gi_id = int(sseqid.split("|")[1])
                gb_acc = sseqid.split("|")[3]
                sseq = sseq.replace("-", "") #TODO here is where we want to grab the full sequence
                sscinames = sscinames.replace(" ", "_").replace("/", "_")
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                stitle = salltitles
                # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
                # Get all of them! (get redundant seq info)
                found_taxids = set()
                found_spn = set()
                if len(sallseqid.split(";")) > 1:
                    # get additional info only for seq that pass the eval
                    if evalue < float(self.config.e_value_thresh):
                        staxids_l = staxids.split(";")
                        sscinames_l = sscinames.split(";")
                        sallseqid_l = sallseqid.split(";")
                        salltitles_l = salltitles.split("<>")
                        count = 0
                        spn_range = 0
                        stop_while = False
                        while len(found_taxids) < len(staxids_l):  # as long as i have not found all taxids for the seq
                            count += 1
                            if stop_while:
                                break
                            if count == 5:
                                break  # too many tries to find correct number of redundant taxa
                            elif count == 1:
                                for i in range(0, len(sallseqid_l)):
                                    if len(found_taxids) == len(staxids_l):
                                        break
                                    gi_id = sallseqid_l[i].split("|")[1]
                                    gb_acc = sallseqid_l[i].split("|")[3]
                                    # if gb acc was already read in before stop the for loop
                                    if gb_acc in query_dict or gb_acc in self.data.gb_dict:
                                        stop_while = True
                                        break
                                    stitle = salltitles_l[i]
                                    # spn_title are used to figure out if gb_acc is from same tax_id
                                    # if both var are the same, we do not need to search GB for taxon info
                                    spn_title_before = salltitles_l[i-1].split(" ")[0:spn_range]
                                    spn_title = salltitles_l[i].split(" ")[0:spn_range]
                                    # sometimes if multiple seqs are merged,
                                    # we lack information about which taxon is which gb_acc...
                                    # test it here:
                                    # if we have same number of ids and taxon id go ahead as usual
                                    #TODO if we grab the full sequence using the accession number, we can also associate it with the correct taxon info 
                                    if len(sallseqid_l) == len(staxids_l):
                                        staxids = staxids_l[i]
                                        sscinames = sscinames_l[i]
                                    # only one taxon id present, all are from same taxon
                                    elif len(staxids_l) == 1:
                                        staxids = staxids_l[0]
                                        sscinames = sscinames_l[0]
                                    elif i != 0 and spn_title != spn_title_before:
                                        read_handle = self.ids.entrez_efetch(gb_acc)
                                        sscinames =  ncbi_data_parser.get_ncbi_tax_name(read_handle).replace(" ", "_").replace("/", "_")
                                        staxids =  ncbi_data_parser.get_ncbi_tax_id(read_handle)
                                        spn_range = len(sscinames.split("_"))
                                    elif i == 0:  # for first item in redundant data, always get info
                                        read_handle = self.ids.entrez_efetch(gb_acc)
                                        sscinames =  ncbi_data_parser.get_ncbi_tax_name(read_handle).replace(" ", "_").replace("/", "_")
                                        staxids = ncbi_data_parser.get_ncbi_tax_id(read_handle)
                                        spn_range = len(sscinames.split("_"))
                                    else:  # if spn_titles were the same, we do not need to add same seq again
                                        continue
                                    assert str(staxids) in staxids_l, (staxids, staxids_l)
                                    # next vars are used to stop loop if all taxids were found
                                    found_taxids.add(staxids)
                                    found_spn.add(sscinames)
                                    if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                                        query_dict[gb_acc] = \
                                            {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                             'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                             'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                            # same loop as above, only that it mostly used entrez searches for tax_id
                            # this is as sometimes the if above does not yield in stop_while == True,
                            # through different taxa names
                            elif count >= 1 and stop_while is False:
                                for i in range(0, len(sallseqid_l)):
                                    if len(found_taxids) == len(staxids_l):
                                        break
                                    gi_id = sallseqid_l[i].split("|")[1]
                                    gb_acc = sallseqid_l[i].split("|")[3]
                                    stitle = salltitles_l[i]
                                    # if gb acc was already read in before stop the for loop
                                    if gb_acc in query_dict or gb_acc in self.data.gb_dict:
                                        stop_while = True
                                        break
                                    # sometimes if multiple seqs are merged, we lack the information
                                    # which taxon is which gb_acc...test it here:
                                    # if we have same number information go ahead as usual
                                    if len(sallseqid_l) == len(staxids_l):
                                        staxids = staxids_l[i]
                                        sscinames = sscinames_l[i]
                                    elif len(staxids_l) == 1:  # only one taxon id present, all are from same taxon
                                        staxids = staxids_l[0]
                                        sscinames = sscinames_l[0]
                                    else:
                                        read_handle = self.ids.entrez_efetch(gb_acc)
                                        sscinames = ncbi_data_parser.get_ncbi_tax_name(read_handle).replace(" ", "_").replace("/", "_")
                                        staxids = ncbi_data_parser.get_ncbi_tax_id(read_handle)
                                        spn_range = len(sscinames.split("_"))
                                    assert str(staxids) in staxids_l, (staxids, staxids_l)
                                    # next vars are used to stop loop if all taxids were found
                                    found_taxids.add(staxids)
                                    found_spn.add(sscinames)
                                    if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                                        query_dict[gb_acc] = \
                                             {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                              'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                              'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                else:
                    staxids = int(staxids)
                    self.ids.spn_to_ncbiid[sscinames] = staxids
                    if gb_acc not in self.ids.acc_ncbi_dict:  # fill up dict with more information.
                        self.ids.acc_ncbi_dict[gb_acc] = staxids
                    if gb_acc not in query_dict and gb_acc not in self.newseqs_acc:
                        query_dict[gb_acc] = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                              'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                              'bitscore': bitscore, 'sseq': sseq, 'title': stitle}
                
        for key in query_dict.keys():
            if float(query_dict[key]["evalue"]) < float(self.config.e_value_thresh):
                gb_acc = query_dict[key]["accession"]
                if len(gb_acc.split(".")) >= 2:
                    if gb_acc not in self.data.gb_dict:
                        self.new_seqs[gb_acc] = query_dict[key]["sseq"]
                        self.data.gb_dict[gb_acc] = query_dict[key]
                # else:
                    # debug("was added before")
            else:

                fn = open("{}/blast_threshold_not_passed.csv".format(self.workdir), "a+")
                fn.write("blast_threshold_not_passed:\n")
                fn.write("{}, {}, {}\n".format(query_dict[key]["sscinames"], query_dict[key]["accession"],
                                               query_dict[key]["evalue"]))
                fn.close()

    def read_unpublished_blast_query(self):
        """
        Reads in the blast files generated during local_blast_for_unpublished() and adds seq to self.data.gb_dict and
        self.new_seqs.

        """
        debug("read unpublished blast query")
        output_blast = "output_tst_fn.xml"
        gb_counter = 1
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        # with cd(os.path.join(self.workdir, "blast")):
        xml_file = open(output_blast)
        os.chdir(general_wd)
        blast_out = NCBIXML.parse(xml_file)
        fn = open("{}/not_added_local_seq.csv".format(self.workdir), "a")
        fn.write("not_added_local_seq")
        for blast_record in blast_out:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    local_id = alignment.title.split("|")[-1].split(" ")[-1]
                    if float(hsp.expect) < float(self.config.e_value_thresh):
                        if local_id not in self.data.gb_dict:  # skip ones we already have
                            unpbl_local_id = "unpubl_{}".format(local_id)
                            self.new_seqs[unpbl_local_id] = hsp.sbjct
                            # debug(self.new_seqs[unpbl_local_id])
                            self.data.gb_dict[unpbl_local_id] = {'title': "unpublished", 'localID': local_id}
                            # debug(self.data.unpubl_otu_json)
                            # debug(local_id)
                            # debug(type(local_id))
                            # debug('otu{}'.format(local_id.replace("_", "").replace("-", "")))
                            self.data.gb_dict[unpbl_local_id].update(
                                self.data.unpubl_otu_json['otu{}'.format(local_id.replace("_", "").replace("-", ""))])
                            gb_counter += 1
                            # debug(self.data.gb_dict[unpbl_local_id])
                            # debug(some)
                    else:
                        fn.write("{}: {}".format(alignment.title.split("|")[-1].split(" ")[-1], hsp.expect))
                        # if local_id not in self.gb_not_added:
                        #     self.gb_not_added.append(local_id)
                        writeinfofiles.write_not_added_info(self, local_id, "threshold not passed")
                        # needs to be deleted from gb_dict,
                        # maybe we find a better fitting blast query seq and then it might get added
                        del self.data.gb_dict[unpbl_local_id]  # print(some)
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from unpublished database\n".format(len(self.new_seqs)))

    def read_webbased_blast_query(self, fn_path):
        """ Implementation to read in results of web blast searches.

        :param fn_path: path to file containing the local blast searches
        :return: updated self.new_seqs and self.data.gb_dict dictionaries
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
                            gb_id = alignment.title.split("|")[3]  # 1 is for gi
                            if len(gb_id.split(".")) == 1:
                                debug(gb_id)
                            if gb_id not in self.data.gb_dict:  # skip ones we already have
                                # gb_id = int(alignment.title.split('|')[1])  # 1 is for gi
                                # assert type(gb_id) is int
                                # SHOULD NOT BE NECESSARY....IS WEBBLAST HAS THE TAXON ALREADY LIMITED
                                # if len(self.acc_list_mrca) >= 1 and (gb_id not in self.acc_list_mrca):
                                #     pass
                                # else:
                                taxid,taxname, seq = self.ids.get_tax_seq_acc(gb_id)
                                self.new_seqs[gb_id] = seq
                                gi_id = alignment.title.split('|')[1]
                                gb_acc = alignment.__dict__['accession']
                                stitle = alignment.__dict__['title']
                                hsps = alignment.__dict__['hsps']
                                length = alignment.__dict__['length']
                                query_dict = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'title': stitle,
                                              'length': length, 'hsps': hsps}
                                self.data.gb_dict[gb_id] = query_dict
                        else:
                            # if gb_id not in self.gb_not_added:
                            #     self.gb_not_added.append(gb_id)
                            #     writeinfofiles.write_not_added_info(self, gb_id, "threshold not passed")
                            writeinfofiles.write_not_added_info(self, gb_id, "threshold not passed")
                            # needs to be deleted from gb_dict,
                            # maybe we find a better fitting blast query seq and then it might get added
                            del self.data.gb_dict[gb_id]
        except ValueError:
            sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

    def read_blast_wrapper(self, blast_dir=None):
        """reads in and processes the blast xml files

        :param blast_dir: path to directory which contains blast files
        :return: fills different dictionaries with information from blast files
        """
        debug("read_blast_wrapper")
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
            self.read_unpublished_blast_query()
        else:
            if not self._blasted:
                self.run_blast_wrapper()
            assert os.path.exists(self.blast_subdir)
            for taxon in self.data.aln:
                # debug(self.config.blast_loc)
                fn = None
                if self.config.blast_loc == "local":
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                if self.config.gb_id_filename is True: #TODO what is this doing?
                    fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label) 
                    if fn is None:
                        fn = self.data.otu_dict[taxon.label].get('^user:TaxonName', taxon.label)
                    fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
                else:
                    fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                if _DEBUG:
                    sys.stdout.write("attempting to read {}\n".format(fn_path))
                if os.path.isfile(fn_path):
                    if self.config.blast_loc == 'local':  # new method to read in txt format
                        self.read_local_blast_query(fn_path)
                    else:
                        self.read_webbased_blast_query(fn_path)
        self.date = str(datetime.date.today())
#        debug("len new seqs dict after evalue filter")
#        debug(len(self.new_seqs))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from GenBank after evalue filtering\n".format(len(self.new_seqs)))

        self._blast_read = 1

    def get_sp_id_of_otulabel(self, label): 
        #TODO This was doing a bunch of searches that sould have already happened when creating the OTU.
        #If we didn't find an ncbi_id then, we shouldn't look again.

        """Get the species name and the corresponding ncbi id of the otu.

        :param label: otu_label = key from otu_dict
        :return: ncbi id of corresponding label
        """
        # debug("get_tax_id_of_otulabel")
        ncbi_id = self.data.otu_dict[label].get("^ncbi:taxon",0)
#        debug("otu {} has taxon id {}".format(label, ncbi_id))
        return ncbi_id


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
        if label == None: #in case of add_otu failure, doean't edit dict 
            debug("otu_id None")
            sys.stderr.write("otu_id None")
            return seq_dict
        id_of_label = self.get_sp_id_of_otulabel(label)
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        continue_search = False
        never_add = False
        for tax_lab in tax_list:
            existing_id = self.get_sp_id_of_otulabel(tax_lab)
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(new_seq) >= sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * self.config.maxlen:
                debug("seq not added because it's too long...")
                self.data.otu_dict[label]['^physcraper:status'] = "not added; sequence is too long"
                gb_id = self.data.otu_dict[label]["^ncbi:accession"]
                # if gb_id not in self.gb_not_added:
                #     self.gb_not_added.append(gb_id)
                #TODO merge with other label search
                if id_of_label in self.ids.ncbiid_to_spn.keys():
                    tax_name = self.ids.ncbiid_to_spn[id_of_label]
                else:
                    tax_name = self.ids.ncbi_parser.get_name_from_id(id_of_label)
                cutoff = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen) * self.config.maxlen
                reason = "sequence too long: new seq len ({}) vs.  cutoff ({})".format(len(new_seq), cutoff)
                writeinfofiles.write_not_added(id_of_label, tax_name, gb_id, reason, self.workdir)
                # writeinfofiles.write_not_added_info(self, local_id, "threshold not passed")
                # needs to be deleted from gb_dict,
                # maybe we find a better fitting blast query seq and then it might get added
#                del self.data.gb_dict[gb_id]
                                    
            elif len(inc_seq) >= len(new_seq):  # if seq is identical and shorter
                if inc_seq.find(new_seq) != -1:
                    if type(existing_id) == int and existing_id != int(id_of_label):  # different otus, add
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, "
                                             "but different species name\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; " \
                                                                          "subsequence, but different taxon"
                        seq_dict[label] = seq
                        debug("{} and {} are subsequences, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    else:  # subseq of same otu
                        if _VERBOSE:
                            sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "subsequence, not added"
                        debug("{} not added, subseq of {}".format(id_of_label, existing_id))
                        never_add = True
                        continue
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
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of {}, but different "
                                             "taxon\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added; supersequence, " \
                                                                          "but different taxon"
                        seq_dict[label] = seq
                        debug("{} and  {} supersequence, but different sp. concept".format(id_of_label, existing_id))
                        continue_search = True
                        continue
                    elif self.data.otu_dict[tax_lab].get('^physcraper:status') == "local seq":
                        if _VERBOSE:
                            sys.stdout.write("seq {} is supersequence of local seq {}, "
                                             "both kept in alignment\n".format(label, tax_lab))
                        self.data.otu_dict[label]['^physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        debug("{} and  {} added".format(id_of_label, existing_id))
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
                        self.data.otu_dict[tax_lab]['^physcraper:status'] = "deleted, {} is supersequence".format(label)
                        debug("{} added, instead of  {}".format(id_of_label, existing_id))
                        self.del_superseq.add(tax_lab)
                        continue_search = True
                        continue

        if continue_search is True or never_add is True:
            if (self.data.otu_dict[label]['^physcraper:status'].split(' ')[0] in self.seq_filter) or never_add is True:
                if label in seq_dict.keys():
                    del seq_dict[label]
                if label in self.data.aln.taxon_namespace or label in self.data.tre.taxon_namespace:
                    self.data.remove_taxa_aln_tre(label)
                # should not be necessary, information what happened to seq should have been added in lines before
                # else:
                #     debug("label was never added to aln or tre")
                # Note: should not be the word 'deleted', as this is used in self.seq_filter
                # self.data.otu_dict[label]['^physcraper:status'] = "removed in seq dict build"
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
        # need to re-calculate orig_seq_len before using it
        self.data.orig_seqlen = [len(self.data.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in
                                 self.data.aln]
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        assert self.config.seq_len_perc <= 1, \
            ("your config seq_len_param is not smaller than 1: {}".format(self.config.seq_len_perc))
        seq_len_cutoff = avg_seqlen * self.config.seq_len_perc
        self.del_superseq = set()  # will contain deleted superseqs for the assert below
        # all_added_gi is to limit the adding to new seqs, if we change the rank of filtering later
        all_added_gi = set()
        for gb_id, seq in self.new_seqs.items():
            if seq == None:
                sys.stderr.write("No sequence found for accession number {}\n".format(gb_id))
                continue
            added = False
            assert gb_id in self.data.gb_dict.keys(), (gb_id, self.data.gb_dict.keys())
            if gb_id not in all_added_gi:
                all_added_gi.add(gb_id)
                # debug(gb_id)
                if len(gb_id.split(".")) == 1:
                    debug("problem gb_id {}".format(gb_id))
                if self.blacklist is not None and gb_id in self.blacklist:
                    debug("gb_id in blacklist, not added")
                    pass
                elif gb_id in self.newseqs_acc:  # added to increase speed. often seq was found in another blast file
                    debug("passed, was already added")
                    pass
                else:
                    # debug(len(seq.replace("-", "").replace("N", "")))
                    # debug(seq_len_cutoff)
                    if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                        # debug(self.config.blast_loc)
                        if self.config.blast_loc == "local":
                            #in local blast need to check if new seqs are in taxon on interest
                            ncbi_id, tax_name = ncbi_data_parser.get_tax_info_from_acc(gb_id, self.data, self.ids)
                            if self.ids.ncbi_parser.match_id_to_mrca(ncbi_id, self.mrca_ncbi):  # belongs to ingroup mrca -> add to data
                                self.newseqs_acc.append(gb_id)
                                otu_id = self.data.add_otu(gb_id, self.ids)
                                self.seq_dict_build(seq, otu_id, tmp_dict)
                                added = True
                                    # debug(some)
                            else:
                                debug("does not match")
                                reason = "not_part_of_mrca: {}".format(self.mrca_ncbi)
                        else: #in remote we have already checked for MRCA
                            added = True
                            self.newseqs_acc.append(gb_id)
                            otu_id = self.data.add_otu(gb_id, self.ids)
                            self.seq_dict_build(seq, otu_id, tmp_dict)
                        
                    else:
                        #debug("seq too short")
                        # do not add them to not added, as seq len depends on sequence which was blasted and
                        # there might be a better matching seq (one which will return a longer sequence...)
                        # if gb_id not in self.gb_not_added:
                        #     self.gb_not_added.append(gb_id)
                        # writeinfofiles.write_not_added_info(self, gb_id, "seqlen_threshold_not_passed")
                        len_seq = len(seq.replace("-", "").replace("N", ""))
                        reason = "seqlen_threshold_not_passed: len seq: {}, min len: {}".format(len_seq, seq_len_cutoff)
 

                # sequences that could not be added are written out to file
                if added is False:
                    if "sscinames" in self.data.gb_dict[gb_id]:
                            tax_name = self.data.gb_dict[gb_id]['sscinames']
                            ncbi_id = self.data.gb_dict[gb_id]['staxids']
                    else:
                        tax_name = None
                        ncbi_id = None
                    writeinfofiles.write_not_added(ncbi_id, tax_name, gb_id, reason, self.workdir)

                        # fn = open("{}/not_added_seq.csv".format(self.workdir), "a+")
                        # fn.write(
                        #     "seqlen_threshold_not_passed, {}, {}, min len: {}\n".format(gb_id, len_seq, seq_len_cutoff))
                        # fn.close()
        # this assert got more complicated, as sometimes superseqs are already deleted in seq_dict_build() ->
        # then subset assert it not True
        old_seqs_ids = set()
        for tax in old_seqs:
            old_seqs_ids.add(tax)
        tmp_dict_plus_super = set()
        for item in tmp_dict.keys():
            tmp_dict_plus_super.add(item)
        for item in self.del_superseq:
            tmp_dict_plus_super.add(item)
        # print(tmp_dict_plus_super)
        # assert old_seqs_ids.issubset(tmp_dict.keys()), ([x for x in old_seqs_ids if x not in tmp_dict.keys()])
        assert old_seqs_ids.issubset(tmp_dict_plus_super), ([x for x in old_seqs_ids if x not in tmp_dict_plus_super])
        for tax in old_seqs:
            if tax in tmp_dict:
                del tmp_dict[tax]
        self.new_seqs_otu_id = tmp_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        debug("len new seqs dict after remove identical")
        debug(len(self.new_seqs_otu_id))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from Genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def dump(self, filename=None, recursion=100000):
        """writes out class to pickle file.
        We need to increase the recursion depth here, as it currently fails with the standard run.

        :param filename: optional filename
        :param recursion: pickle often failed with recursion depth, that's why it's increased here
        :return: writes out file
        """
        current = sys.getrecursionlimit()
        sys.setrecursionlimit(recursion)

        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("{}/scrape_checkpoint.p".format(self.workdir), "wb")
        pickle.dump(self, ofi, pickle.HIGHEST_PROTOCOL)
        sys.setrecursionlimit(current)

    def write_query_seqs(self):
        """writes out the query sequence file"""
        debug("write query seq")
        if not self._blast_read:
            self.read_blast_wrapper()
        self.newseqs_file = "{}.fasta".format(self.date)
        fi = open("{}/{}".format(self.workdir, self.newseqs_file), "w")
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
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        if _VERBOSE:
            sys.stdout.write("aligning query sequences \n")
        self.data._reconcile()  # I think reconcile is what was needed here...instead of alien hack
        # note: sometimes there are still sp in any of the aln/tre
        # hack for the alien taxa thing
        self.remove_alien_aln_tre()
        self.data.write_papara_files()
        os.chdir(self.workdir)  # Clean up dir moving
        # with cd(self.workdir):
        assert self.data.aln.taxon_namespace == self.data.tre.taxon_namespace
        try:
            subprocess.check_call(["papara",
                                   "-t", "random_resolve.tre",
                                   "-s", "aln_ott.phy",
                                   #  "-j", "{}".format(self.config.num_threads),  # FIXME: Does not work on some machines
                                   "-q", self.newseqs_file,
                                   "-n", papara_runname])  # FIXME directory ugliness
            if _VERBOSE:
                sys.stdout.write("Papara done")
        except subprocess.CalledProcessError as grepexc:
            print "error code", grepexc.returncode, grepexc.output
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it installed?\n")
                sys.exit(-5)
            # handle file not found error.
            else:
                # Something else went wrong while trying to run `wget`
                raise
        path = "{}/papara_alignment.{}".format(self.workdir, papara_runname)
        assert os.path.exists(path), "{path} does not exists".format(path=path)
        os.chdir(cwd)
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment."
                                                    "{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        self._query_seqs_aligned = 1

    def remove_alien_aln_tre(self):
        """Sometimes there were alien entries in self.tre and self.aln.

        This function ensures they are properly removed."""

        treed_tax = set()
        for leaf in self.data.tre.leaf_nodes():
            treed_tax.add(leaf.taxon)
        aln_tax = set()
        for tax, seq in self.data.aln.items():
            aln_tax.add(tax)
        prune = treed_tax ^ aln_tax
        del_tre = []
        del_aln = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_tax)
            if taxon in aln_tax:
                # debug(taxon)
                del_aln.append(taxon)
            if taxon in treed_tax:
                del_tre.append(taxon)
        # debug(del_aln)
        # debug(del_tre)
        self.data.aln.remove_sequences(del_aln)
        self.data.tre.prune_taxa(del_tre)
        for tax_lab in self.data.aln.taxon_namespace:
            if tax_lab not in self.data.tre.taxon_namespace:
                sys.stderr.write("tax {} not in tre. This is an alien name in the data.\n".format(tax_lab))
                self.data.remove_taxa_aln_tre(tax_lab)
        for tax_lab in self.data.tre.taxon_namespace:
            if tax_lab not in self.data.aln.taxon_namespace:
                sys.stderr.write("tax {} not in aln. This is an alien name in the data.\n".format(tax_lab))
                self.data.remove_taxa_aln_tre(tax_lab)
        # this should not need to happen here
        # self.data.prune_short()
        # self.data.trim()

    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree."""
        if self.backbone is True:
            with cd(self.workdir):
                backbonetre = Tree.get(path="{}/backbone.tre".format(self.workdir),
                                       schema="newick",
                                       preserve_underscores=True)

                backbonetre.resolve_polytomies()
                backbonetre.write(path="random_resolve.tre", schema="newick", unquoted_underscores=True)

        if os.path.exists("RAxML_labelledTree.PLACE"):
            os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
        if _VERBOSE:
            sys.stdout.write("placing query sequences \n")
        with cd(self.workdir):
            try:
                debug("try")
                subprocess.call(["raxmlHPC-PTHREADS",
                                 "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-f", "v",
                                 "-s", "papara_alignment.extended",
                                 "-t", "random_resolve.tre",
                                 "-n", "PLACE"])
                placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                    schema="newick",
                                    preserve_underscores=True)
            except:
                try:
                    subprocess.call(["raxmlHPC",
                                     "-m", "GTRCAT",
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
                        sys.exit(-6)
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
        self._query_seqs_placed = 1

    def est_full_tree(self, path="."):
        """Full raxml run from the placement tree as starting tree.
        The PTHREAD version is the faster one, hopefully people install it. if not it falls back to the normal raxml.
        the backbone options allows to fix the sceleton of the starting tree and just newly estimates the other parts.
        """
        debug("est full tree")
        cwd = os.getcwd()
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[-1]))
        try:
            num_threads = int(self.config.num_threads)
            if self.backbone is not True:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
            else:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")

            if self.backbone is not True:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
            else:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(self.date)])
        os.chdir(cwd)
        self._full_tree_est = 1

    def calculate_bootstrap(self):
        """Calculates bootstrap and consensus trees.

        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        -z: specifies file with multiple trees

        """
        debug("calculate bootstrap")
        cwd = os.getcwd()
        os.chdir(self.workdir)
        # with cd(self.workdir):
        # # check if job was started with mpi
        # # this checks if the whole file was started as mpiexec
        # env_var = [os.environ.get('PMI_RANK'), os.environ.get('PMI_SIZE'), os.environ.get('OMPI_COMM_WORLD_SIZE')]
        # mpi = False
        # for var in env_var:
        #     if var is not None:
        #         mpi = True
        # check if job was started with mpi
        # this checks if actual several cores and nodes were allocated
        ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
        nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
        # env_var = int(nnodes) * int(ntasks)
        mpi = False
        if nnodes is not None and ntasks is not None:
            env_var = int(nnodes) * int(ntasks)
            mpi = True
        if mpi:
            debug("run with mpi")
            subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2",
                             # "raxmlHPC-PTHREADS", "-T", "{}".format(num_threads),
                             "-m", "GTRCAT",
                             "-s", "previous_run/papara_alignment.extended",
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "{}".format(self.date)])
        else:
            try:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-s", "previous_run/papara_alignment.extended",
                                 "-p", "1", "-b", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])
            except:
                subprocess.call(["raxmlHPC",
                                 "-m", "GTRCAT",
                                 "-s", "previous_run/papara_alignment.extended",
                                 "-p", "1", "-b", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])

        try:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-s", "previous_run/papara_alignment.extended",
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "all{}".format(self.date)])

            # strict consensus:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "STRICT",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "StrictCon{}".format(self.date)])
            # majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MR",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "MR_{}".format(self.date)])
            # extended majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MRE",
                             "-z", "RAxML_bootstrap.all{}".format(self.date),
                             "-n", "EMR{}".format(self.date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")
            # make bipartition tree
            # is the -f b command
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
        os.chdir(cwd)

    def remove_blacklistitem(self):
        """This removes items from aln, and tree, if the corresponding Genbank identifer were added to the blacklist.

        Note, that seq that were not added because they were similar to the one being removed here, are lost
        (that should not be a major issue though, as in a new blast_run, new seqs from the taxon can be added.)
        """
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            acc = self.data.otu_dict[tax.label].get("^ncbi:accession")
            if gi_id in self.blacklist or acc in self.blacklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, Genbank identifier is part of blacklist"
        # this should not need to happen here: prune_short; instead...
        self.data.check_tre_in_aln()
        # self.data.prune_short()
        # debug(self.data.tre.as_string(schema='newick'))

    def generate_streamed_alignment(self):
        """runs the key steps and then replaces the tree and alignment with the expanded ones"""
        debug("generate streamed aln")
        if self.blacklist:
            self.remove_blacklistitem()
        debug(len(self.new_seqs))
        debug(len(self.new_seqs_otu_id))
        if len(self.new_seqs) == 0 or len(self.new_seqs_otu_id) == 0:
            if _VERBOSE:
                sys.stdout.write("No new sequences found.\n")
            # self.repeat = 0
            self.calculate_final_tree()
            self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))
        elif len(self.new_seqs) > 0:
            self.data.write_files()  # should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:  # TODO rename to something more intuitive
                self.data.check_tre_in_aln()
                self.write_query_seqs()
                self.align_query_seqs()
                self.place_query_seqs()
                self.data.prune_short()
                self.data.trim()
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
                if self.config.gb_id_filename is not True:
                    os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    if not os.path.exists("{}/previous_run".format(self.workdir)):
                        os.makedirs('{}/previous_run/'.format(self.workdir))
                    if self.config.gb_id_filename is not True:
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                    else:
                        os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[-1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file),
                          "{}/previous_run/newseqs.fasta".format(self.workdir))
                self.data.write_labelled(label='^physcraper:TaxonName', add_gb_id=True)
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {}  # Wipe for next run
                self.new_seqs_otu_id = {}
                self.repeat = 1
            else:
                self.calculate_final_tree()
                self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))
        #     else:
        #         if _VERBOSE:
        #             sys.stdout.write("No new sequences after filtering.\n")
        #         # self.repeat = 0
        #         self.calculate_final_tree()
        #         self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))

        # else:
        #     if _VERBOSE:
        #         sys.stdout.write("No new sequences found.\n")
        #     # self.repeat = 0
        #     self.calculate_final_tree()
        #     self.data.dump("{}/final_ATT_checkpoint.p".format(self.workdir))

        self.reset_markers()

        filter_by_local_blast.del_blastfiles(self.workdir)  # delete local blast db
        self.data.dump()
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))

    def calculate_final_tree(self):
        """Calculates the final tree using a trimmed alignment.

        :return: final PS data
        """
        debug("calculate final tree")
        self.data.write_files(treepath="physcraper_final_notrim.tre", alnpath="physcraper_final_notrim.fas")
        self.data.prune_short()
        self.data.trim()
        self.data.write_files(treepath="physcraper_final_trim.tre", alnpath="physcraper_final_trim.fas")
        if os.path.exists("[]/previous_run".format(self.workdir)):
            self.est_full_tree(path="previous_run")
        else:
            self.est_full_tree()
        self.repeat = 0
        self.calculate_bootstrap()

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
                filter_by_local_blast.write_filterblast_db(self.workdir, key, seq, fn="local_unpubl_seq")
        with cd(os.path.join(self.workdir, "blast")):
            cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
            os.system(cmd1)

