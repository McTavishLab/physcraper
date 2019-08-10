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

from copy import deepcopy

from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from physcraper.configobj import ConfigObj
from physcraper.ids import IdDicts
from physcraper.aligntreetax import AlignTreeTax
from physcraper.helpers import cd, get_raxml_ex, write_filterblast_db
from physcraper.ncbi_data_parser import get_gi_from_blast, get_acc_from_blast

from physcraper import ncbi_data_parser
from physcraper import writeinfofiles
from physcraper import AWSWWW

_VERBOSE = 1
_DEBUG = 0
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
    def __init__(self, data_obj, ids_obj, ingroup_mrca=None, threshold = 5):
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
        self.write_mrca()
        self.threshold = threshold
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

    def write_mrca(self):
        with open('{}/mrca.txt'.format(self.workdir),"w") as fi:
            fi.write("ingroup mrca ott_id {}\n".format(self.mrca_ott))
            fi.write("ingroup mrca NCBI: {}\n".format(self.mrca_ncbi))


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
        assert os.path.isdir(self.config.blastdb), ("blast dir does not exist: '{}'.".format(self.config.blastdb))
        with cd(self.config.blastdb):
            # this format (6) allows to get the taxonomic information at the same time
            outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid'"
            # outfmt = " -outfmt 5"  # format for xml file type
            # TODO query via stdin
            blastcmd = "blastn -query " + "{}/tmp.fas".format(abs_blastdir) + \
                       " -db {}/nt -out ".format(self.config.blastdb) + abs_fn + \
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
            debug("blasting {} using {}".format(fn_path, self.config.url_base))
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
        today = str(datetime.date.today()).replace("-", "/")
        debug("run_blast_wrapper")
        debug(self.blast_subdir)
        self._blast_read = 0
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
#            tmpfile = open("{}/{}.tmp".format(self.workdir, otu_id),"w")
#            tmpfile.write(seq.symbols_as_string())
#            tmpfile.write("\n")
            assert otu_id in self.data.otu_dict
            last_blast = self.data.otu_dict[otu_id].get('^physcraper:last_blasted')
            if last_blast == None:
                time_passed = delay + 1
            else:
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(
                last_blast, "%Y/%m/%d")).days)
            if time_passed > delay:
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
               # tmpfile.write(query)
                if self.config.blast_loc == "local":
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
                # if _DEBUG:
                #     sys.stdout.write("attempting to write {}\n".format(fn_path))
                if not os.path.isfile(fn_path):
                    if _VERBOSE:
                        sys.stdout.write("blasting seq {}\n".format(taxon.label))
                    if self.config.blast_loc == 'local':
                        self.run_local_blast_cmd(query, taxon.label, fn_path)
                    if self.config.blast_loc == 'remote':
                        if last_blast:
                            equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi, last_blast, today)
                        else:
                            equery = "txid{}[orgn]".format(self.mrca_ncbi)
                        debug(equery)
               #         tmpfile.write("\nequery\n")
               #         tmpfile.write(equery)
               #         tmpfile.write("\nquery\n")
               #         tmpfile.write(query)
               #         tmpfile.close()
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
#        debug("read_local_blast_query")
        query_dict = {}
        with open(fn_path, mode="r") as infile:
            for lin in infile:
                sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, stitle, sallseqid = lin.strip().split('\t')
                gb_acc = get_acc_from_blast(sseqid)
                gi_id = get_gi_from_blast(sseqid)                
                sseq = sseq.replace("-", "") #TODO here is where we want to grab the full sequence MK: I wrote a batch query for the seqs we are interested. Makes it faster.
                taxname = sscinames.replace(" ", "_").replace("/", "_")
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                if gb_acc.split('.')>= 2: 
                    # get additional info only for seq that pass the eval
                    if evalue < float(self.config.e_value_thresh):
                        if gb_acc not in self.new_seqs.keys(): # do not do it for gb_ids we already considered
                            # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
                            # Get all of them! (they can be of a different taxon ids = get redundant seq info)
                            if len(sallseqid.split(';')) > 1:
                                debug("multiple taxa have identical blast match {}, using {}\n".format(sallseqid, gb_acc))
                                taxid, taxname, full_seq = self.ids.get_tax_seq_acc(gb_acc)
                            else:
                                taxid = int(staxids)
                                self.ids.spn_to_ncbiid[sscinames] = staxids
                                if gb_acc not in self.ids.acc_ncbi_dict:  # fill up dict with more information.
                                    self.ids.acc_ncbi_dict[gb_acc] = staxids
                                full_seq = self.get_full_seq(gb_acc, sseq)
                            query_dict[gb_acc] = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'staxids': staxids,
                                                          'sscinames': sscinames, 'pident': pident, 'evalue': evalue,
                                                          'bitscore': bitscore, 'sseq': full_seq, 'title': stitle}
                            self.data.gb_dict[gb_acc] = query_dict[gb_acc]
                            self.new_seqs[gb_acc] = query_dict[gb_acc]["sseq"]

                    else:
                        fn = open("{}/blast_threshold_not_passed.csv".format(self.workdir), "a+")
                        fn.write("blast_threshold_not_passed: {}, {}, {}\n".format(sscinames, gb_acc, gi_id))
                        fn.close()
                else:
                    pass
                    #sys.stdout.write("skipping acc {}, unexpected format\n".format(gb_acc))



 

    def get_full_seq(self, gb_acc, blast_seq, local=False):
        """
        Get full sequence from gb_acc that was retrieved via blast. 

        Currently only used for local searches, Genbank database sequences are retrieving them in batch mode, which is hopefully faster.

        :param gb_acc: unique sequence identifier (often genbank accession number)
        :param blast_seq: sequence retrived by blast,
        :return: full sequence, the whole submitted sequence, not only the part that matched the blast query sequence
        """
        # debug("get full seq")
        # if we did not already try to get full seq:

        if not os.path.exists("{}/tmp".format(self.workdir)):
            os.mkdir("{}/tmp".format(self.workdir))
        if not os.path.exists("{}/tmp/full_seq_{}".format(self.workdir, gb_acc)):
            fn = "{}/tmp/tmp_search.csv".format(self.workdir)
            fn_open = open(fn, "w+")
            fn_open.write("{}\n".format(gb_acc.split(".")[0]))
            fn_open.close()
            db_path = "{}/nt".format(self.config.blastdb)

            cmd1 = "blastdbcmd -db {}  -entry_batch {} -outfmt %f -out {}/tmp/full_seq_{}.fasta".format(db_path, fn, self.workdir, gb_acc)
            # debug(cmd1)
            os.system(cmd1)
            # read in file to get full seq        
            fn = "{}/tmp/full_seq_{}.fasta".format(self.workdir, gb_acc)
            f = open(fn)
            seq = ""
            for line in iter(f):
                line = line.rstrip().lstrip()
                if line[0]  != ">":
                    seq += line
                elif line[0]  == ">":
                    pass
#                    assert gb_acc in line, (gb_acc, line)
            f.close()
        # check direction of sequence:
        orig = Seq(seq, generic_dna)
        dna_comp = orig.complement()
        dna_rcomp = orig.reverse_complement()
        dna_r = orig[::-1]
        full_seq = str()
        if blast_seq.replace("-", "") in orig:
            full_seq = seq
        elif blast_seq.replace("-", "") in dna_r:
            full_seq = dna_r
        elif blast_seq.replace("-", "") in dna_comp:
            full_seq = dna_comp
        elif blast_seq.replace("-", "") in dna_rcomp:
            full_seq = dna_rcomp
        if blast_seq.replace("-", "") not in full_seq:
            taxid, taxname, full_seq = self.ids.get_tax_seq_acc(gb_acc)
        full_seq = str(full_seq)
        assert type(full_seq) == str, (type(full_seq))
        return full_seq


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
                                 sys.stdout.write("skipping acc {}, incorrect format\n".format(gb_id))
                            elif gb_id not in self.data.gb_dict:  # skip ones we already have
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
                            writeinfofiles.write_not_added_info(self, gb_id, "evalue threshold not passed")
                            # needs to be deleted from gb_dict,
                            # maybe we find a better fitting blast query seq and then it might get added
                            #del self.data.gb_dict[gb_id]
        except ValueError:
            sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

    def read_blast_wrapper(self, blast_dir=None):
        """reads in and processes the blast xml files

        :param blast_dir: path to directory which contains blast files
        :return: fills different dictionaries with information from blast files
        """
#        debug("read_blast_wrapper")
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
                self.run_blast_wrapper()
        assert os.path.exists(self.blast_subdir)
        for taxon in self.data.aln:
                # debug(self.config.blast_loc)
            if self.config.blast_loc == "local":
                file_ending = "txt"
            else:
                file_ending = "xml"
#                if self.config.gb_id_filename is True: #TODO what is this doing?
#                    fn = self.data.otu_dict[taxon.label].get('^ncbi:accession', taxon.label) 
#                    if fn is None:
#                        fn = self.data.otu_dict[taxon.label].get('^user:TaxonName', taxon.label)
#                    fn_path = "{}/{}.{}".format(self.blast_subdir, fn, file_ending)
#                else:
            fn_path = "{}/{}.{}".format(self.blast_subdir, taxon.label, file_ending)
            if _DEBUG:
                sys.stdout.write("reading {}\n".format(fn_path))
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
        if len(self.new_seqs) == 0:
            sys.stderr.write("no new squences found in blast. Exiting")
            sys.exit()
        self.remove_identical_seqs()
        self._blast_read = 1


    def seq_dict_build(self, seq, new_otu_label, seq_dict):
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
        #debug("new_lab: {} in seq_dict_build".format(new_otu_label))
        if new_otu_label == None: #in case of add_otu failure, doean't edit dict 
            return seq_dict
        tax_new_seq = self.data.otu_dict[new_otu_label].get('^ncbi:taxon', 1)
        if self.config.blast_loc == "local": #need to check if included in taxon of intrest (mrca)
            if self.ids.ncbi_parser.match_id_to_mrca(tax_new_seq, self.mrca_ncbi):
                #taxon is within mrca, continue
                #debug("local: otu {} is within mrca".format(new_otu_label))
                pass
            else:
                #debug("local: otu {} is NOT within mrca".format(new_otu_label))
                self.data.otu_dict[new_otu_label]['^physcraper:status'] = "Not within requesdted MRCA ncbi:{}".format(self.mrca_ncbi)
                return seq_dict
        #debug("otu {} has tax_id {}".format(new_otu_label, tax_new_seq))
        new_seq = seq.replace("-", "").lower()
        otu_list = deepcopy(seq_dict.keys())
        should_add = True
        reason = 'new'
        i = 0
        assert new_otu_label not in otu_list
        for otu_lab in otu_list:
            #debug("old lab: {}".format(otu_lab))
            i += 1
            if _VERBOSE:
                sys.stdout.write(".")
                if i % 50 == 0:
                    sys.stdout.write("\n")
            existing_tax_id = self.data.otu_dict[otu_lab].get('^ncbi:taxon', None)
            inc_seq = seq_dict[otu_lab].replace("-", "").lower()
            if len(inc_seq) >= len(new_seq):  
                #debug("seq {} is shorter than {}".format(new_otu_label, otu_lab))
                if new_seq in inc_seq:# if seq is identical and shorter
                    #debug("seq is identical and shorter")
                    if existing_tax_id != tax_new_seq:  # different taxa
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write("seq {} is subsequence of {}, "
                                             "but different species name ".format(new_otu_label, otu_lab))
                        reason = "new seq added; subsequence of {}, but different taxon".format(otu_lab)
                        #debug("otu {}, tmp status {}".format(new_otu_label, reason))
                        #still should be added, but need to check other samples
                    else:  # subseq of same otu
                        reason = "seq {} is subsequence or identical to {}, not added ".format(new_otu_label, otu_lab)
                        if _VERBOSE:
                            sys.stdout.write(reason)
                        self.data.otu_dict[new_otu_label]['^physcraper:status'] = "subsequence, not added"
                        #debug("{} not added, subseq of {}".format(new_otu_label, otu_lab))
                        #debug("{} was NOT added to seq_dict: {}".format(new_otu_label, reason))
                        return seq_dict
                else:
                    pass
                    #didn't run into any problems yet, should_add still true
            elif len(new_seq) > len(inc_seq):  
                #debug("seq is longer")
                if new_seq.find(inc_seq) != -1:
                    if self.data.otu_dict[otu_lab].get('^physcraper:status') == "original":
                        reason = "seq {} is supersequence of original seq {}, "\
                                             "both kept in alignment ".format(new_otu_label, otu_lab)
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write("\n" + reason)
                    elif existing_tax_id != tax_new_seq:  # different taxa
                        reason = "seq {} is supersequence of {}, but different taxon ".format(new_otu_label, otu_lab)
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write("\n" + reason)
                        #can still be added
                    else:
                        # new seq s a super sequence, delet old one and add new one. DO NOT KEEP CHECKING
                        del seq_dict[otu_lab]
                        seq_dict[new_otu_label] = seq
                        self.data.remove_taxa_aln_tre(otu_lab)
                        reason = "\nseq {} is supersequence of {}, {} added and {} removed ".format(new_otu_label, otu_lab, new_otu_label, otu_lab)
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write(reason)
                        self.data.otu_dict[otu_lab]['^physcraper:status'] = "deleted, {} is supersequence ".format(new_otu_label)
                        self.data.otu_dict[new_otu_label]['^physcraper:status'] = "new seq added in place of {} ".format(otu_lab)
                        seq_dict[new_otu_label] = seq
                        self.data.otu_dict[new_otu_label]['^physcraper:status'] = reason
                        #debug("{} was added to seq_dict: {}".format(new_otu_label, reason))
                        return seq_dict
        seq_dict[new_otu_label] = seq
        self.data.otu_dict[new_otu_label]['^physcraper:status'] = reason
        #debug("{} was added to seq_dict: {}".format(new_otu_label, reason))
        return seq_dict

    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        """
        #debug("remove identical seqs")
        #debug("new seqs keys are {}".format(self.new_seqs.keys()))
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
        seq_len_min = avg_seqlen * self.config.seq_len_perc
        seq_len_max = avg_seqlen * self.config.maxlen
        #debug("we already have {}".format(all_added_gi))
        for gb_id, seq in self.new_seqs.items():
            assert gb_id in self.data.gb_dict.keys()
#            debug("gb_id is {}".format(gb_id) )
            assert seq
            if seq_len_min < len(seq) < seq_len_max:
                if self.blacklist is not None and gb_id in self.blacklist:
                    debug("gb_id {} in blacklist, not added".format(gb_id))
                    pass
                else:
                    otu_id = self.data.add_otu(gb_id, self.ids)
                    tmp_dict = self.seq_dict_build(seq, otu_id, tmp_dict)
            else:
                debug("len {}:{} was not between {} and {}".format(gb_id, len(seq), seq_len_min, seq_len_max))
        otu_in_aln = set([taxon.label for taxon in self.data.aln])
        for otu in otu_in_aln:
            del tmp_dict[otu]
        filter_dict = self.filter_seqs(tmp_dict, type='random', threshold = self.threshold)
        self.new_seqs_otu_id = filter_dict  # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        self.new_seqs = {} #Wipe clean
#        debug("len new seqs otu dict after remove identical{}".format(len(self.new_seqs_otu_id)))
        sys.stdout.write("**** Found {} new sequences****\n".format(len(self.new_seqs_otu_id)))
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from Genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()
    
    def filter_seqs(self, tmp_dict, threshold=5, type="random"):
        assert type in ['length','random'], "type {} not recognized, please filter by 'length' or 'random'".format(type)
        selected_otus = set()
        filtered_dict = {}
        new_sp_d = self.make_sp_dict(tmp_dict.keys())
        debug("There are {} taxa in the new taxa".format(len(new_sp_d)))
        debug("The keys of tmp_dict".format(tmp_dict.keys()))
        aln_otus = set([taxon.label for taxon in self.data.aln])
        aln_sp_d = self.make_sp_dict(aln_otus)
        debug("There are {} taxa in aln".format(len(aln_sp_d)))
        alltax = set(new_sp_d.keys()).union(aln_sp_d.keys())
        sys.stdout.write("{} taxa in orginal alignment; {} taxa in updated alignemnt {}, keeping max {} seq per taxon".format(len(aln_sp_d), len(alltax)), threshold)
        for tax_id in new_sp_d:
            debug(" {} new seqs for taxon {}".format(len(new_sp_d[tax_id]), tax_id))
            tax_otus = []
            current_num = len(aln_sp_d.get(tax_id, []))
            if current_num > threshold:
#                sys.stdout.write("no sequences added for taxon {}, as {} already in alignmenet".format(tax_id, current_num))
                pass #already enough
            else:
                count = threshold - len(aln_sp_d.get(tax_id,[]))
                otu_list = new_sp_d[tax_id]
                if count > len(otu_list):
                    tax_otus = otu_list
                else:
                    if type == 'random':
                        tax_otus = self.select_seq_at_random(otu_list, count)
                    if type == 'length':
                        tax_otus = self.select_seq_by_length(otu_list, tmp_dict, count)
                debug("passing on {} otus for tax id {}".format(len(tax_otus), tax_id))
                assert isinstance(selected_otus, set), 'why not set?!?!'
                assert isinstance(tax_otus, list), 'why not list?!?!'
                selected_otus.update(tax_otus)
        for otu in tmp_dict:
            if otu in selected_otus:
                filtered_dict[otu] = tmp_dict[otu]
        return filtered_dict



    def make_sp_dict(self, otu_list=[], downtorank=None):
        """Mkaes dict of OT_ids by species"""
        self.downtorank = downtorank
        if otu_list == []:
            otu_list = self.new_seqs_otu_id.keys()
        debug("make sp_dict")
        sp_d = {}
        for otu_id in otu_list:
            if self.data.otu_dict[otu_id]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                tax_id = self.data.otu_dict[otu_id].get('^ncbi:taxon') 
                if tax_id == None:
                    tax_id = "X"
                if self.downtorank is not None:
                    downtorank_name = None
                    downtorank_id = None
                    if self.config.blast_loc == 'remote':
                        sys.stderr.write("Filtering by taxon ranks not functional for remote ncbi searches yet.")
                        sys.exit(-7)
                        #tax_id = self.ids.get_rank_info_from_web(taxon_name=tax_name)
                        #lineage2ranks = self.ids.otu_rank[tax_id]["rank"]
                        #ncbi = NCBITaxa()
                        #if lineage2ranks == 'unassigned':
                        #    downtorank_id = tax_id
                        #    downtorank_name = tax_name
                        #else:
                        #    for key_rank, val in lineage2ranks.items():
                        #        if val == downtorank:
                        #            downtorank_id = key_rank
                        #            value_d = ncbi.get_taxid_translator([downtorank_id])
                        #            downtorank_name = value_d[int(downtorank_id)]
                    else:
                        downtorank_id = self.ids.ncbi_parser.get_downtorank_id(tax_id, self.downtorank)
                        downtorank_name = self.ids.ncbi_parser.get_name_from_id(downtorank_id)
                    tax_name = downtorank_name
                    tax_id = downtorank_id
                if tax_id in sp_d:
                    sp_d[tax_id].append(otu_id)
                else:
                    sp_d[tax_id] = [otu_id]
        return sp_d

    def select_seq_by_length(self, otu_list, seq_dict, count):
        """This is another mode to filter the sequences, if there are more than the threshold.

        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".
        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program how many
        sequences for the taxon are already available in the aln.

        !!! sometimes the only seq in seq_w_maxlen is the original seq,
        then this is the one to be added, but it will be removed,
        later as it is no new seq! thus no new seq for that species is added

        :param taxon_id: key from self.sp_seq_d
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        debug("select_seq_by_length")
        lens = []
        otu_len_dict = {}
        for otu in otu_list:
            otu_len_dict[otu] = len(seq_dict[otu])
            lens.append(len(seq_dict[otu]))
        lens.sort(reverse=True)
        cutoff = lens[count]
#        debug("cutoff is {}".format(cutoff))
#        debug("lengths is {}".format(lens))
        selected_otus = []
        for otu in otu_len_dict:
#            debug("otu {} has len {}".format(otu, otu_len_dict[otu]))
            if otu_len_dict[otu] >= cutoff:
#                debug("selected!")
                selected_otus.append(otu)
                if len(selected_otus) == count:
                    return selected_otus
        return selected_otus
        

    def select_seq_at_random(self, otu_list, count):
        """This is another mode to filter the sequences, if there are more than the threshold.

        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".
        Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program how many
        sequences for the taxon are already available in the aln.

        !!! sometimes the only seq in seq_w_maxlen is the original seq,
        then this is the one to be added, but it will be removed,
        later as it is no new seq! thus no new seq for that species is added

        :param taxon_id: key from self.sp_seq_d
        :param count: self.count_num_seq(tax_id)["seq_present"]
        :return: self.filtered_seq
        """
        debug("select_seq_at random")
        sample_count = min(count, len(otu_list))
        selected_otus = random.sample(otu_list, sample_count)
        return selected_otus



    def pickle_dump(self, filename=None, recursion=100000):
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

    def dump(self):
        self.data.write_otus('otu_dict.json')
        with open("{}/{}".format(self.workdir, 'config.json'), "w") as outfile:
                json.dump(self.config.__dict__, outfile)

    def write_query_seqs(self, filename='date'):
        """writes out the query sequence file"""
        debug("write query seq")
        if not self._blast_read:
            self.read_blast_wrapper()
        if filename == 'date':
            self.newseqs_file = "{}.fasta".format(self.date)
        else:
            self.newseqs_file = filename
        fipath = ("{}/{}".format(self.workdir, self.newseqs_file))
        if _VERBOSE:
            sys.stdout.write("writing out sequences\n")
        with open(fipath, "w") as fi:
            for otu_id in self.new_seqs_otu_id.keys():
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1
        return fipath

    def write_all_unaligned(self, filename='date'):
        """writes out the query sequence file"""
        debug("write query + aligned seqs")
        if self._blasted == 1 and self._blast_read == 0:
            self.read_blast_wrapper()
        if filename == 'date':
            finam = "{}_ALL.fasta".format(self.date)
        else:
            finam = filename
        fipath =  "{}/{}".format(self.workdir, finam)  
        #write out existing
        self.data.aln.write(path=fipath, schema='fasta')
        #append new
        if _VERBOSE:
            sys.stdout.write("writing out ALL sequences\n")
        with open(fipath, "a") as fi:
            for otu_id in self.new_seqs_otu_id.keys():
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        fi.close()
        self._query_seqs_written = 1
        self.data.write_otus(schema='table')
        return fipath


    def align_query_seqs(self, aligner = 'muscle'):
        if not self._blast_read:
            self.read_blast_wrapper()
        assert aligner in ['muscle', 'papara']
        if aligner == 'papara':
            self.run_papara()
        if aligner == 'muscle':
            self.run_muscle()
        self.data.trim()
        alnfi = self.data.write_aln()
        return alnfi
    
    def replace_aln(self, filename, schema = 'fasta'):
        newaln = DnaCharacterMatrix.get(path=filename, schema=schema)
        for taxon in newaln:
            assert taxon.label in self.data.otu_dict
        self.data.aln = newaln
        self.new_seqs_otu_id = {}
        self._blasted = 0
        self.blast_read = 0


    def replace_tre(self, filename, schema = 'newick'):
        newtre= Tree.get(path=filename,
                   schema=schema,
                   taxon_namespace = self.data.aln.taxon_namespace)
        aln_tax = [taxon.label for taxon in self.data.aln]
        for taxon in newtre.leaf_nodes():
            assert taxon.taxon.label in self.data.otu_dict
            assert taxon.taxon.label in aln_tax
        self.data.tre = newtre



    def run_muscle(self, outname = 'muscle_aln.fas'):
        finame = self.write_all_unaligned()
        outpath = "{}/{}".format(self.workdir, outname)
        try:
            subprocess.check_call(["muscle",
                                   "-in", finame,
                                   "-out", outpath])
            if _VERBOSE:
                sys.stdout.write("Muscle done")
        except subprocess.CalledProcessError as grepexc:
            print "error code", grepexc.returncode, grepexc.output
        self.replace_aln(outpath)


    def run_papara(self, papara_runname="extended"):
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
        aln = DnaCharacterMatrix.get(path="{}/papara_alignment."
                                                    "{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = True  # Was too strict...
        if _VERBOSE:
            sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        return aln
        self._query_seqs_aligned = 1



    def place_query_seqs(self, alignment = None, query = None, tree = "random_resolve.tre"):
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
        rac_ex = get_raxml_ex()
        with cd(self.workdir):
            try:
                debug("try")
                subprocess.call([rax_ex,
                                 "-m", "GTRCAT",
                                 "-f", "v",
                                 "-s", alignment,
                                 "-t", tree,
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
                        raise
            placetre.resolve_polytomies()
            placetre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True)
        self._query_seqs_placed = 1


    def est_full_tree(self, alignment = None, startingtree = None, backbone = False, method = "raxml"):
        """Full raxml run from the placement tree as starting tree.
        The PTHREAD version is the faster one, hopefully people install it. if not it falls back to the normal raxml.
        the backbone options allows to fix the sceleton of the starting tree and just newly estimates the other parts.
        """
        cwd = os.getcwd()
        if alignment == None:
            debug("call align query seqs from est full tree, self._blast_read is {}".format(self._blast_read))
            alignment = self.align_query_seqs()
        if startingtree == None:
            startingtree = os.path.abspath(self.data.write_random_resolve_tre())
        debug("est full tree")
        os.chdir(self.workdir)
        rax_ex = get_raxml_ex()
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/treest_prev".format(self.workdir))
        num_threads = int(self.config.num_threads)
        label = "{}".format(self.date)
        if self.backbone:
            cmd= [rax_ex, "-m", "GTRCAT", "-s", alignment, "-r", "backbone.tre", "-p", "1", "-n", label]
            
        else:
            cmd= [rax_ex, "-m", "GTRCAT", "-s", alignment, "-t", startingtree, "-p", "1", "-n", label]
        process = subprocess.Popen(cmd)
        process.wait()
        if _VERBOSE:
                sys.stdout.write("running: "+" ".join(cmd)+"\n")
        self.replace_tre("RAxML_bestTree.{}".format(label))
        self.data.write_labelled('^ot:ottTaxonName')
        os.chdir(cwd)
        return label

#        self.tre = Tree.get(data="RAxML_bestTree.{}".format(label),
#                            schema="newick",
#                            preserve_underscores=True,
#                            taxon_namespace=self.data.aln.taxon_namespace)
#        self.data._reconcile()
#        self.data.write_files()


                    
    def calculate_bootstrap(self, alignment = None):
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
        if alignment == None:
            alignment = "previous_run/papara_alignment.extended"
        with cd(self.workdir):
            ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
            nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
            # env_var = int(nnodes) * int(ntasks)
            mpi = False
            if nnodes is not None and ntasks is not None:
                env_var = int(nnodes) * int(ntasks)
                mpi = True
            rax_ex = get_raxml_ex()
            if mpi:
                debug("run with mpi")
                subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2",
                                 # "raxmlHPC-PTHREADS", "-T", "{}".format(num_threads),
                                 "-m", "GTRCAT",
                                 "-s", alignment,
                                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])
            else:
                subprocess.call([rax_ex, "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-s", alignment,
                                 "-p", "1", "-b", "1", "-#", "autoMRE",
                                 "-n", "{}".format(self.date)])

            # try:
            #     subprocess.call([rax_ex, "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
            #                      "-s", alignment,
            #                      "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
            #                      "-n", "all{}".format(self.date)])

            #     # strict consensus:
            #     subprocess.call([rax_ex,  "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
            #                      "-J", "STRICT",
            #                      "-z", "RAxML_bootstrap.all{}".format(self.date),
            #                      "-n", "StrictCon{}".format(self.date)])
            #     # majority rule:
            #     subprocess.call([rax_ex, "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
            #                      "-J", "MR",
            #                      "-z", "RAxML_bootstrap.all{}".format(self.date),
            #                      "-n", "MR_{}".format(self.date)])
            #     # extended majority rule:
            #     subprocess.call([rax_ex, "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
            #                      "-J", "MRE",
            #                      "-z", "RAxML_bootstrap.all{}".format(self.date),
            #                      "-n", "EMR{}".format(self.date)])

            

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
                write_filterblast_db(self.workdir, key, seq, fn="local_unpubl_seq")
        with cd(os.path.join(self.workdir, "blast")):
            cmd1 = "makeblastdb -in {}_db -dbtype nucl".format("local_unpubl_seq")
            os.system(cmd1)

