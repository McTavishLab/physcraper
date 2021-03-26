"""
The core blasting and new sequence integration module
"""
import sys
import os
import subprocess
import datetime
import glob
import random
from dendropy import Tree, DnaCharacterMatrix
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from physcraper.ids import IdDicts
from physcraper.aligntreetax import AlignTreeTax
from physcraper.helpers import cd, get_raxml_ex
from physcraper.ncbi_data_parser import get_gi_from_blast, get_acc_from_blast
from physcraper.opentree_helpers import root_tree_from_synth
from physcraper import AWSWWW

_VERBOSE = 0

def set_verbose():
    """Set output to verbose"""
    global _VERBOSE # pylint: disable=global-statement
    _VERBOSE = 1

_DEBUG = 0
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)


class PhyscraperScrape():
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
          * **self.new_seqs_otu_id**: dictionary that contains
          the new sequences that passed the remove_identical_seq() step:

            * key: otu_id
            * value: see otu_dict, is a subset of the otu_dict, all sequences that will be newly added to aln and tre
          * **self.mrca_ncbi**:  int ncbi identifier of mrca

          * **self.blast_subdir**: path to folder that contains the files writen during blast

          * **self.newseqs_file**: filename of files that contains the sequences from self.new_seqs_otu_id
          * **self.date**: Date of the run - may lag behind real date!
          * **self.repeat**: either 1 or 0, it is used to determine if we continue updating the tree,
          no new seqs found = 0
          * **self.newseqs_acc**: list of all gi_ids that were passed into remove_identical_seq().
          Used to speed up adding process
          * **self.blocklist**: list of gi_id of sequences that shall not be added or need to be removed.
          Supplied by user.
          * **self.seq_filter**: list of words that may occur in otu_dict.status and which shall not be
          used in the building of FilterBlast.sp_d
           (that's the main function), but it is also used as assert statement to make sure unwanted seqs are not added.
          * **self.unpublished**: True/False. Used to look for local unpublished seq that shall be added if True.
          * **self.path_to_local_seq:** Usually False, contains path to unpublished sequences if option is used.

        Following functions are called during the init-process:

            * **self.reset_markers()**: adds things to self: I think they are used to make sure certain function run,
            if program crashed and pickle file is read in.
                * self._blasted: 0/1, if run_blast_wrapper() was called, it is set to 1 for the round.
                * self._blast_read: 0/1, if read_blast_wrapper() was called, it is set to 1 for the round.
                * self._identical_removed: 0
                * self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
                * self._query_seqs_aligned: 0
                * self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
                * self._reconciled: 0
                * self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
    """
    def __init__(self, data_obj, ids_obj=None, search_taxon=None):
        assert isinstance(data_obj, AlignTreeTax)
        self.data = data_obj
        self.workdir = data_obj.workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.inputsdir = "{}/inputs_{}".format(self.workdir, self.data.tag)
        if not os.path.exists(self.inputsdir):
            os.makedirs(self.inputsdir)
        self.data.write_files(direc=self.inputsdir)
        self.outputsdir = "{}/outputs_{}".format(self.workdir, self.data.tag)
        if not os.path.exists(self.outputsdir):
            os.makedirs(self.outputsdir)
        self.rundir = "{}/run_{}".format(self.workdir, self.data.tag)
        if not os.path.exists(self.rundir):
            os.makedirs(self.rundir)
        self.logfile = "{}/logfile".format(self.rundir)
        if ids_obj is None:
            self.ids = IdDicts()
        else:
            assert isinstance(ids_obj, IdDicts)
            self.ids = ids_obj
        self.config = self.ids.config  # pointer to config
        self.config.write_file(self.rundir)
        self.new_seqs = {}  # all new seq after read_blast_wrapper
        self.new_seqs_otu_id = {}  # only new seq which passed remove_identical
        self.blast_subdir = "{}/blast_run_{}".format(self.workdir, self.data.tag)

        self.date = str(datetime.date.today()) #Date of the run-may lag behind real date!
        self.repeat = 1  # used to determine if we continue updating the tree
        self.newseqs_acc = []  # all ever added Genbank accession numbers during any PhyScraper run
        self.seq_filter = ['deleted', 'subsequence,', 'not', "removed", "deleted,",
                           "local"]
        self.reset_markers()
        self.gb_not_added = []  # list of blast seqs not added
        self.del_superseq = set()  # items that were deleted bc they are superseqs, needed for assert statement
        self.mrca_ott = data_obj.mrca_ott
        self.mrca_ncbi = self.ids.ott_to_ncbi.get(data_obj.mrca_ott)
        if self.mrca_ncbi is None or self.mrca_ncbi == 0:
            sys.stderr.write("search mrca{} does not have a direct match to ncbi.\n".format(search_taxon))
        sys.stdout.write("restricting blast search\
                         to taxon ncbi:{} (ott:{}; {})\n".format(self.mrca_ncbi,
                                                                 self.mrca_ott,
                                                                 self.ids.ott_to_name[self.mrca_ott]))
        debug("created physcraper ncbi_mrca {},".format(self.mrca_ncbi))
        self.map_taxa_to_ncbi()
        assert self.mrca_ncbi
        self.write_mrca()
        self.data.write_labelled(filename="taxonname", label='^ot:ottTaxonName', direc=self.inputsdir)
        self.data.write_otus(schema='table', direc=self.inputsdir)
        if not os.path.exists("{}/otu_info_{}.json".format(self.rundir, self.data.tag)):
            self.data.write_otus(schema='json', direc=self.rundir)
        self.threshold = self.config.spp_threshold
#markers for status
#holders for new data
        self.blocklist = []

    def map_taxa_to_ncbi(self):
        """Find NCBI ids for taxa from OpenTree"""
        for otu in self.data.otu_dict:
            if self.data.otu_dict[otu].get("^ncbi:taxon") is None:
                if self.data.otu_dict[otu].get("^ot:ottId"):
                    ottid = self.data.otu_dict[otu]["^ot:ottId"]
                    self.data.otu_dict[otu]["^ncbi:taxon"] = self.ids.ott_to_ncbi.get(ottid, 0)


    def reset_markers(self):
        """set completion markers back to 0 for a re-run"""
        self._blasted = 0
        self._blast_read = 0
        self._query_seqs_written = 0
        self._query_seqs_placed = 0
        self._full_tree_est = 0

    def write_mrca(self):
        """Write out search info to file"""
        with open('{}/mrca.txt'.format(self.inputsdir), "w") as fi:
            fi.write("search mrca ott_id {}\n".format(self.mrca_ott))
            fi.write("search mrca NCBI: {}\n".format(self.mrca_ncbi))


    def run_local_blast_cmd(self, query, taxon_label, fn_path):
        """Contains the cmds used to run a local blast query, which is different from the web-queries.

        :param query: query sequence
        :param taxon_label: corresponding taxon name for query sequence
        :param fn_path: path to output file for blast query result

        :return: runs local blast query and writes it to file
        """
        abs_blastdir = os.path.abspath(self.blast_subdir)
        abs_outfile = os.path.abspath(fn_path)
        toblast = open("{}/tmp.fas".format(os.path.abspath(self.blast_subdir)), "w+")
        toblast.write(">{}\n".format(taxon_label))
        toblast.write("{}\n".format(query))
        toblast.close()
        assert os.path.isdir(self.config.blastdb), ("blast dir does not exist: '{}'.".format(self.config.blastdb))
        # this format (6) allows to get the taxonomic information at the same time
        outfmt = "6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid"
        # outfmt = " -outfmt 5"  # format for xml file type
        blastcmd = ["blastn",
                    "-query",
                    "{}/tmp.fas".format(abs_blastdir),
                    "-db",
                    "{}/nt".format(self.config.blastdb),
                    "-out",
                    abs_outfile,
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    str(self.config.num_threads),
                    "-max_target_seqs",
                    str(self.config.hitlist_size),
                    "-max_hsps",
                    str(self.config.hitlist_size)]
        try:
            subprocess.check_output(blastcmd, cwd=self.config.blastdb)
        except KeyboardInterrupt:
            if os.stat(abs_outfile).st_size == 0:
                os.remove(abs_outfile)
            sys.stderr.write("KeyboardInterrupt")
            sys.exit()



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

    def run_blast_wrapper(self):
        """generates the blast queries and saves them depending on the blasting method to different file formats

        It runs blast if the sequences was not blasted since the user defined threshold in the config file (delay).

        :return: writes blast queries to file
        """
        delay = self.config.delay
        today = str(datetime.date.today()).replace("-", "/")
        debug("Today's date is {}".format(today))
        debug("run_blast_wrapper")
        debug(self.blast_subdir)
        debug("current alignment length {}".format(len(self.data.aln)))
        self._blast_read = 0
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        #try:
        for taxon, seq in self.data.aln.items():
            sys.stdout.write("Blasting {}\n".format(taxon))
            otu_id = taxon.label
            assert otu_id in self.data.otu_dict
            last_blast = self.data.otu_dict[otu_id].get('^physcraper:last_blasted')
            if last_blast is None:
                time_passed = delay + 1
            else:
                time_passed = abs((datetime.datetime.strptime(today, "%Y/%m/%d") -
                                   datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days)
            # following line should stay == False, because sometimes ingroup status is just unknown
            if self.data.otu_dict[otu_id].get("^physcraper:ingroup") == False: # pylint: disable=singleton-comparison
                sys.stdout.write("tip {} not in ingroup. Will not blast, \n".format(otu_id))
                continue
            if time_passed > delay:
                query = seq.symbols_as_string().replace("-", "").replace("?", "")
                if self.config.blast_loc == "local":
                    file_ending = "txt"
                else:
                    file_ending = "xml"
                idtag = "{}{}".format(taxon.label, self.data.otu_dict[taxon.label].get('^ncbi:accession', ''))
                fn_path = "{}/{}.{}".format(self.blast_subdir, idtag, file_ending)
                # if _DEBUG:
                #     sys.stdout.write("attempting to write {}\n".format(fn_path))
                if not os.path.isfile(fn_path):
                    if self.config.blast_loc == 'local':
                        self.run_local_blast_cmd(query, taxon.label, fn_path)
                    if self.config.blast_loc == 'remote':
                        if last_blast:
                            equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi, last_blast, today)
                        else:
                            equery = "txid{}[orgn]".format(self.mrca_ncbi)
                        debug(equery)
                        self.run_web_blast_query(query, equery, fn_path)
                    self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                else:
                    t = os.path.getmtime(fn_path)
                    filedate = datetime.date.fromtimestamp(t)
                    self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = str(filedate).replace("-", "/")
                    if _DEBUG:
                        sys.stdout.write("file {} exists in current blast run. Will not blast, "
                                         "delete file to force\n".format(fn_path))
            else:
                sys.stdout.write("otu {} was last blasted on {}, {} days ago and is not being re-blasted."
                                     "Use run_blast_wrapper(delay = 0) to force a search.\n".format(otu_id,
                                                                                                    last_blast,
                                                                                                    time_passed))
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
                if gb_acc is None:
                    continue
                gi_id = get_gi_from_blast(sseqid)
                sseq = sseq.replace("-", "")
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                if len(gb_acc.split('.')) > 1:
                    # get additional info only for seq that pass the eval
                    if evalue < float(self.config.e_value_thresh):
                        if gb_acc not in self.new_seqs.keys(): # do not do it for gb_ids we already considered
                            # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
                            # Get all of them! (they can be of a different taxon ids = get redundant seq info)
                            if len(sallseqid.split(';')) > 1:
                                for match in sallseqid.split(';'):
                                    gb_acc = get_acc_from_blast(match)
                                    if gb_acc is not None:
                                        full_seq = self.get_full_seq(gb_acc, sseq)
                                        query_dict[gb_acc] = {'^ncbi:gi': gi_id,
                                                              'accession': gb_acc,
                                                              'staxids': staxids,
                                                              'sscinames': sscinames,
                                                              'pident': pident,
                                                              'evalue': evalue,
                                                              'bitscore': bitscore,
                                                              'sseq':full_seq,
                                                              'title': stitle}
                                        self.data.gb_dict[gb_acc] = query_dict[gb_acc]
                                        self.new_seqs[gb_acc] = query_dict[gb_acc]["sseq"]
                            else:
                                self.ids.spn_to_ncbiid[sscinames] = staxids
                                if gb_acc not in self.ids.acc_ncbi_dict:  # fill up dict with more information.
                                    self.ids.acc_ncbi_dict[gb_acc] = staxids
                                full_seq = self.get_full_seq(gb_acc, sseq)
                                query_dict[gb_acc] = {'^ncbi:gi': gi_id,
                                                      'accession': gb_acc,
                                                      'staxids': staxids,
                                                      'sscinames': sscinames,
                                                      'pident': pident,
                                                      'evalue': evalue,
                                                      'bitscore': bitscore,
                                                      'sseq': full_seq,
                                                      'title': stitle}
                                self.data.gb_dict[gb_acc] = query_dict[gb_acc]
                                self.new_seqs[gb_acc] = query_dict[gb_acc]["sseq"]

                    else:
                        fn = open("{}/blast_threshold_not_passed.csv".format(self.rundir), "a+")
                        fn.write("blast_threshold_not_passed: {}, {}, {}, {}\n".format(sscinames,
                                                                                       gb_acc,
                                                                                       gi_id,
                                                                                       evalue))
                        fn.close()
                else:
                    pass
                    #sys.stdout.write("skipping acc {}, unexpected format\n".format(gb_acc))





    def get_full_seq(self, gb_id, blast_seq):
        """
        Get full sequence from gb_acc that was retrieved via blast.

        Currently only used for local searches, Genbank database sequences are retrieving them in batch mode,
        which is hopefully faster.

        :param gb_acc: unique sequence identifier (often genbank accession number)
        :param blast_seq: sequence retrived by blast,
        :return: full sequence, the whole submitted sequence, not only the part that matched the blast query sequence
        """
        # debug("get full seq")
        # if we did not already try to get full seq:

        if not os.path.exists(self.ids.full_seq_path):
            os.mkdir(self.ids.full_seq_path)
        seq_path = "{}/{}.fasta".format(self.ids.full_seq_path, gb_id)
        if not os.path.exists(seq_path):
            db_path = "{}/nt".format(self.config.blastdb)
            try:
                subprocess.check_call(["blastdbcmd",
                                       "-db", db_path,
                                       "-entry", gb_id,
                                       "-outfmt", "%f",
                                       "-out", seq_path])

            except subprocess.CalledProcessError as grepexc:
                sys.stderr.write("error code {}, {},\n".format(grepexc.returncode, grepexc.output))
                sys.stderr.write("downloading from NCBI\n")
                read_handle = Entrez.efetch(db="nucleotide", id=gb_id, retmode="xml")
                seq = read_handle[0][u'GBSeq_sequence']
                with open(seq_path, 'w') as fi:
                    fi.write("> {} \n".format(gb_id))
                    fi.write(seq)
                fi.close()
                #sys.exit()
            # read in file to get full seq
        f = open(seq_path)
        seq = ""
        for line in iter(f):
            line = line.rstrip().lstrip()
            if line[0] != ">":
                seq += line
            elif line[0] == ">":
                pass
#                    assert gb_acc in line, (gb_acc, line)
        f.close()
        # check direction of sequence:
        match = blast_seq.replace("-", "")
        full_seq = self.check_complement(match, seq, gb_id)
        #if blast_seq.replace("-", "") not in full_seq:
        #    taxid, taxname, full_seq = self.ids.get_tax_seq_acc(gb_acc)
        full_seq = str(full_seq)
        assert isinstance(full_seq, str), (type(full_seq))
        return full_seq

    def check_complement(self, match, seq, gb_id):
        """Double check if blast match is to sequence, complement or reverse complement,
        and return correct seq"""
        match = match.lower()
        seq = seq.lower()
        if match in seq.lower():
            return seq
        if isinstance(seq, str):
            orig = Seq(seq)
        elif isinstance(seq, Seq):
            orig = seq
        else:
            sys.stderr.write("check complement requires a Bio.Seq object\n")
        dna_comp = orig.complement()
        dna_rcomp = orig.reverse_complement()
        dna_r = orig[::-1]
        if match in dna_r:
            with open(self.logfile, "a") as log:
                log.write("acc {}, reverse\n".format(gb_id))
            return str(dna_r)
        if match in dna_comp:
            with open(self.logfile, "a") as log:
                log.write("acc {}, comp\n".format(gb_id))
            return str(dna_comp)
        if match in dna_rcomp:
            with open(self.logfile, "a") as log:
                log.write("acc {}, rcomp\n".format(gb_id))
            return str(dna_rcomp)
        return None

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
                        gb_id = alignment.title.split("|")[3]  # 1 is for gi
                        if float(hsp.expect) < float(self.config.e_value_thresh):
                            if len(gb_id.split(".")) == 1:
                                sys.stdout.write("skipping acc {}, incorrect format\n".format(gb_id))
                            elif gb_id not in self.data.gb_dict:  # skip ones we already have
                                seq = self.ids.get_tax_seq_acc(gb_id)
                                gi_id = alignment.title.split('|')[1]
                                gb_acc = alignment.accession
                                stitle = alignment.title
                                match = hsp.sbjct.lower().replace('-', '')
                                if match not in seq.lower():
                                    seq = self.check_complement(match, seq, gb_id)
                                length = alignment.length
                                query_dict = {'^ncbi:gi': gi_id, 'accession': gb_acc, 'title': stitle,
                                              'length': length, 'hsp': hsp}
                                self.data.gb_dict[gb_id] = query_dict
                                self.new_seqs[gb_id] = seq
                            else:
                                pass
                        else:
                            fi = open("{}/blast_threshold_not_passed.csv".format(self.rundir), "a+")
                            fi.write("{}, {}\n".format(gb_id, hsp.expect))
                            fi.close()
        except ValueError:
            sys.stderr.write("Problem reading {}, skipping\n".format(fn_path))

    def read_blast_wrapper(self, blast_dir=None):
        """reads in and processes the blast xml files

        :param blast_dir: path to directory which contains blast files
        :return: fills different dictionaries with information from blast files
        """
        debug("read_blast_wrapper")
        debug("current alignment length {}".format(len(self.data.aln)))
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
        sys.stdout.write("reading blast data\n")
        assert os.path.exists(self.blast_subdir)
        for taxon in self.data.aln:
                # debug(self.config.blast_loc)
            if self.config.blast_loc == "local":
                file_ending = "txt"
            else:
                file_ending = "xml"
            idtag = "{}{}".format(taxon.label, self.data.otu_dict[taxon.label].get('^ncbi:accession', ''))
            fn_path = "{}/{}.{}".format(self.blast_subdir, idtag, file_ending)
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
            sys.stderr.write("no new sequences found in blast. Exiting")
            sys.exit()
        ret = self.remove_identical_seqs()
        self._blast_read = 1
        return ret


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
        if new_otu_label is None: #in case of add_otu failure, doean't edit dict
            return seq_dict
        tax_new_seq = self.data.otu_dict[new_otu_label].get('^ncbi:taxon', 1)
        if self.config.blast_loc == "local": #need to check if included in taxon of intrest (mrca)
            if self.ids.ncbi_parser.match_id_to_mrca(tax_new_seq, self.mrca_ncbi):
                #taxon is within mrca, continue
                #debug("local: otu {} is within mrca".format(new_otu_label))
                pass
            else:
                #debug("local: otu {} is NOT within mrca".format(new_otu_label))
                self.data.otu_dict[new_otu_label]['^physcraper:status'] = "Not within requested MRCA \
                                                                           ncbi:{}".format(self.mrca_ncbi)
                return seq_dict
        #debug("otu {} has tax_id {}".format(new_otu_label, tax_new_seq))
        new_seq = seq.replace("-", "").lower()
        otu_list = seq_dict.keys()
        reason = 'new'
        i = 0
        if new_otu_label in otu_list: ##Whyy would we hit an otu twice?
            return seq_dict
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
                    if int(existing_tax_id) != int(tax_new_seq):  # different taxa
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write("\nseq {} is identical and/or subsequence of {}, "
                                             "but different species name \n".format(new_otu_label, otu_lab))
                        reason = "new seq added; identical and/or subsequence of {}, \
                                  but different taxon".format(otu_lab)
                        #debug("otu {}, tmp status {}".format(new_otu_label, reason))
                        #still should be added, but need to check other samples
                    else:  # subseq of same otu
                        reason = "seq {} is subsequence or identical to {}, not added ".format(new_otu_label, otu_lab)
                        if _VERBOSE:
                            sys.stdout.write("\n{}\n".format(reason))
                        self.data.otu_dict[new_otu_label]['^physcraper:status'] = "identical and/or subsequence \
                                                                                   of {}, not added".format(otu_lab)
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
                            sys.stdout.write("\n{}\n".format(reason))
                        #can still be added
                    else:
                        # new seq s a super sequence, delet old one and add new one. DO NOT KEEP CHECKING
                        del seq_dict[otu_lab]
                        seq_dict[new_otu_label] = seq
                        self.data.remove_taxa_aln_tre(otu_lab)
                        reason = "seq {} is supersequence of {}, {} added \
                                  and {} removed ".format(new_otu_label, otu_lab, new_otu_label, otu_lab)
                        if _VERBOSE or _DEBUG:
                            sys.stdout.write("\n{}\n".format(reason))
                        self.data.otu_dict[otu_lab]['^physcraper:status'] = "deleted, \
                                                                             {} is supersequence".format(new_otu_label)
                        self.data.otu_dict[new_otu_label]['^physcraper:status'] = "new seq added \
                                                                                   in place of {} ".format(otu_lab)
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
            return None
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        # Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        # need to re-calculate orig_seq_len before using it
        self.data.orig_seqlen = [len(self.data.aln[tax].symbols_as_string().replace("-", "").replace("N", ""))
                                 for tax in self.data.aln]
        avg_seqlen = sum(self.data.orig_seqlen) / len(self.data.orig_seqlen)  # HMMMMMMMM
        seq_len_min = avg_seqlen * self.config.minlen
        seq_len_max = avg_seqlen * self.config.maxlen
        #debug("we already have {}".format(all_added_gi))
        for gb_id, seq in self.new_seqs.items():
            assert gb_id in self.data.gb_dict.keys()
#            debug("gb_id is {}".format(gb_id) )
            assert seq
            if seq_len_min < len(seq) < seq_len_max:
                if self.blocklist is not None and gb_id in self.blocklist:
                    debug("gb_id {} in blocklist, not added".format(gb_id))
                else:
                    otu_id = self.data.add_otu(gb_id, self.ids)
                    tmp_dict = self.seq_dict_build(seq, otu_id, tmp_dict)
            else:
                lr = open("{}/seqlen_mismatch.txt".format(self.outputsdir), "a")
                #taxid, taxname, seq = self.ids.get_tax_seq_acc(gb_id)
                lr.write("acc: {}, len: {}\n".format(gb_id, len(seq)))
                lr.close()
                debug("\nlen {}:{} was not between {} and {}\n".format(gb_id, len(seq), seq_len_min, seq_len_max))
        otu_in_aln = {taxon.label for taxon in self.data.aln}
        for otu in otu_in_aln:
            del tmp_dict[otu]
        filter_dict = self.filter_seqs(tmp_dict, selection='random', threshold=self.config.spp_threshold)
        self.new_seqs_otu_id = filter_dict
        # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        self.new_seqs = {} #Wipe clean
#        debug("len new seqs otu dict after remove identical{}".format(len(self.new_seqs_otu_id)))
        sys.stdout.write("**** Found {} new sequences after \
                          removing identical seqs****\n".format(len(self.new_seqs_otu_id)))
        if len(self.new_seqs_otu_id) == 0:
            return None
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from Genbank after removing identical seq, "
                      "of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        return len(self.new_seqs_otu_id)

    def filter_seqs(self, tmp_dict, selection="random", threshold=None):
        """Subselect from sequences to a threshold of number of seqs per species"""
        if threshold is None:
            threshold = int(self.config.spp_threshold)
        msg = "selection {} not recognized, please filter by 'length' or 'random'"
        assert selection in ['random'], msg.format(selection)
        selected_otus = set()
        filtered_dict = {}
        new_sp_d = self.make_sp_dict(tmp_dict.keys())
        debug("There are {} taxa in the new taxa".format(len(new_sp_d)))
        debug("The keys of tmp_dict are {}".format(tmp_dict.keys()))
        aln_otus = {taxon.label for taxon in self.data.aln}
        aln_sp_d = self.make_sp_dict(aln_otus)
        debug("There are {} taxa in aln".format(len(aln_sp_d)))
        alltax = set(new_sp_d.keys()).union(aln_sp_d.keys())
        sys.stdout.write("{} taxa in orginal alignment; \
                          {} taxa in updated alignemnt, \
                          keeping max {} seq per taxon\n".format(len(aln_sp_d), len(alltax), threshold))
        for tax_id in new_sp_d:
            debug(" {} new seqs for taxon {}".format(len(new_sp_d[tax_id]), tax_id))
            tax_otus = []
            current_num = len(aln_sp_d.get(tax_id, []))
            if current_num > threshold:
                pass #already enough
            else:
                count = threshold - len(aln_sp_d.get(tax_id, []))
                otu_list = new_sp_d[tax_id]
                if count > len(otu_list):
                    tax_otus = otu_list
                else:
                    if selection == 'random':
                        tax_otus = self.select_seq_at_random(otu_list, count)
                debug("passing on {} otus for tax id {}".format(len(tax_otus), tax_id))
                assert isinstance(selected_otus, set), 'why not set?!?!'
                assert isinstance(tax_otus, list), 'why not list?!?!'
                selected_otus.update(tax_otus)
        for otu in tmp_dict:
            if otu in selected_otus:
                filtered_dict[otu] = tmp_dict[otu]
            else:
                self.data.otu_dict[otu]['^physcraper:status'] = "removed in sampling down \
                                                                 to {} per spp.".format(self.threshold)
        return filtered_dict



    def make_sp_dict(self, otu_list=None):
        """Makes dict of OT_ids by species"""
        if otu_list is None:
            otu_list = self.new_seqs_otu_id.keys()
        debug("make sp_dict")
        sp_d = {}
        for otu_id in otu_list:
            if self.data.otu_dict[otu_id]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                tax_id = self.data.otu_dict[otu_id].get('^ncbi:taxon')
                if tax_id is None:
                    tax_id = "X"
                if tax_id in sp_d:
                    sp_d[tax_id].append(otu_id)
                else:
                    sp_d[tax_id] = [otu_id]
            else:
                sys.stdout.write("Removing {},{}".format(otu_id, self.data.otu_dict[otu_id]['^physcraper:status']))
        return sp_d


    def select_seq_at_random(self, otu_list, count):
        """elects sequences at random if there are more than the threshold.
        """
        debug("select_seq_at random")
        sample_count = min(count, len(otu_list))
        selected_otus = random.sample(otu_list, sample_count)
        return selected_otus


    def write_new_seqs(self, filename='date'):
        """writes out the query sequence file"""
        debug("write query seq")
        new_seqs = len(self.new_seqs_otu_id.keys())
        if not self._blast_read:
            new_seqs = self.read_blast_wrapper()
        if new_seqs is None:
            return None
        if filename == 'date':
            self.newseqs_file = "NEW{}_{}.fasta".format(self.date, self.data.tag)
        else:
            self.newseqs_file = filename
        fipath = ("{}/{}".format(self.rundir, self.newseqs_file))
        if _VERBOSE:
            sys.stdout.write("writing out sequences\n")
        with open(fipath, "w") as fi:
            for otu_id in self.new_seqs_otu_id.keys():
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1
        return fipath


    def align_new_seqs(self, aligner='muscle'):
        """Align the new sequences against each other"""
        ret = len(self.new_seqs_otu_id)
        if not self._blast_read:
            ret = self.read_blast_wrapper()
        if ret is None:
            return None
        assert aligner in ['muscle', 'papara']
#        if aligner == 'papara':
#            self.run_papara()
        if aligner == 'muscle':
            self.run_muscle()
        alnfi = self.data.write_aln(direc=self.rundir)
        if alnfi is None:
            return None
        self.data.write_otus(schema='table', direc=self.outputsdir)
        self.data.write_otus(schema='json', direc=self.rundir)
        return alnfi

    def replace_aln(self, filename, schema='fasta'):
        """Replace the alignment in the data object with the new alignment"""
        newaln = DnaCharacterMatrix.get(path=filename, schema=schema)
        for taxon in newaln:
            assert taxon.label in self.data.otu_dict, taxon.label
        debug("updating current alignment from file {}, with {} seqs.".format(filename, len(newaln)))
        self.data.aln = newaln
        self.new_seqs_otu_id = {}
        self._blasted = 0
        self.blast_read = 0



    def replace_tre(self, filename, schema='newick'):
        """Replace the tree in the data object with the new tree"""
        newtre = Tree.get(path=filename,
                          schema=schema,
                          preserve_underscores=True,
                          taxon_namespace=self.data.aln.taxon_namespace)
        try:
            newtre = root_tree_from_synth(newtre, self.data.otu_dict)
        except:
            sys.stderr.write("Tree not rooted, root before running conflict analyses.\n")
        self.data.tre = newtre



    def run_muscle(self, input_aln_path=None, new_seqs_path=None, outname='all_align'):
        """Aligns the new sequences and the profile aligns to the exsiting alignment"""
        debug("running muscle\n")
        outpath_ALL = "{}/{}_{}.fas".format(self.rundir, outname, self.data.tag)
        if os.path.exists(outpath_ALL):
            self.replace_aln(outpath_ALL)
            return outpath_ALL
        if input_aln_path is None:
            aln_filename = "original_{}.fas".format(self.data.tag)
            aln_path = "{}/{}".format(self.inputsdir, aln_filename)
            if os.path.exists(aln_path):
                input_aln_path = aln_path
            else:
                input_aln_path = self.data.write_aln(filename=aln_filename, direc=self.rundir)
        else:
            assert os.path.exists(input_aln_path)
        if new_seqs_path is None:
            new_filename = "new_seqs_UNaligned_{}_{}.fas".format(self.date, self.data.tag)
            tmp_new_seqs_path = "{}/{}".format(self.rundir, new_filename)
            if os.path.exists(tmp_new_seqs_path):
                new_seqs_path = tmp_new_seqs_path
            else:
                new_seqs_path = self.write_new_seqs(filename=new_filename)
        else:
            assert os.path.exists(new_seqs_path)
        outpath_NEW = "{}/new_seqs_aligned_{}_{}.fas".format(self.rundir, self.date, self.data.tag)
        f = open('{}/muscle.log'.format(self.rundir), 'a')
        try:
            subprocess.check_call(["muscle",
                                   "-in", new_seqs_path,
                                   "-out", outpath_NEW], stdout=f, stderr=subprocess.STDOUT)
            if _VERBOSE:
                sys.stdout.write("Muscle NEW done.\n")
        except subprocess.CalledProcessError as grepexc:
            sys.stderr.write("error code {}, {}".format(grepexc.returncode, grepexc.output))
        f = open('{}/muscle.log'.format(self.rundir), 'a')
        try:
            cleaned_align_path = "{}/original_cleaned.fas".format(self.rundir)
            cleaned_align_file = open(cleaned_align_path, 'w')
            subprocess.check_call(["sed", "s/?/-/g",
                                   input_aln_path], stdout=cleaned_align_file, stderr=f)
            cleaned_align_file.close()
            subprocess.check_call(["muscle", "-profile",
                                   "-in1", cleaned_align_path,
                                   "-in2", outpath_NEW,
                                   "-out", outpath_ALL], stdout=f, stderr=subprocess.STDOUT)
            if _VERBOSE:
                sys.stdout.write("Muscle ALL done.\n")
        except subprocess.CalledProcessError as grepexc:
            sys.stderr.write("error code {}, {}".format(grepexc.returncode, grepexc.output))
        self.replace_aln(outpath_ALL)
        return outpath_ALL


    def est_full_tree(self, alignment='default', startingtree=None):
        """Full RAxML run from the placement tree as starting tree.
        The PTHREAD version is the faster one, hopefully people install it if not it falls back to the normal RAxML.
        """
        cwd = os.getcwd()
        if alignment == 'default':
            debug("call align query seqs from est full tree, self._blast_read is {}".format(self._blast_read))
            alignment = self.align_new_seqs()
        if alignment is None:
            return None
        if startingtree is None:
            startingtree = os.path.abspath(self.data.write_random_resolve_tre(direc=self.rundir))
        debug("est full tree")
        os.chdir(self.rundir)
        rax_ex = get_raxml_ex()
        for filename in glob.glob('{}/RAxML*'.format(self.rundir)):
            os.rename(filename, "{}/treest_prev".format(self.rundir))
        # num_threads = int(self.config.num_threads)
        label = "{}".format(self.date)
        cmd = [rax_ex, "-m", "GTRCAT", "-s", alignment, "-t", startingtree, "-p", "1", "-n", label]
        process = subprocess.Popen(cmd)
        process.wait()
        if _VERBOSE:
            sys.stdout.write("running: "+" ".join(cmd)+"\n")
        outfilename = "RAxML_bestTree.{}".format(label)
        os.chdir(cwd)
        besttreepath = "{}/{}".format(self.rundir, outfilename)
        return besttreepath

    def calculate_bootstrap(self, alignment='default', num_reps="100"):
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
        if alignment == 'default':
            debug("call align query seqs from est full tree, self._blast_read is {}".format(self._blast_read))
            alignment = self.align_new_seqs()
        if alignment is None:
            sys.stderr.write("no default alignemnt found for bootstrap\n")
            return None
        with cd(self.rundir):
            ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
            nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
            # env_var = int(nnodes) * int(ntasks)
            mpi = False
            if nnodes is not None and ntasks is not None:
                env_var = int(nnodes) * int(ntasks)
                mpi = True
            rax_ex = get_raxml_ex()
            if mpi:
                ## Currently doesn't run
                debug("run with mpi")
                subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2",
                                 # "raxmlHPC-PTHREADS", "-T", "{}".format(num_threads),
                                 "-m", "GTRCAT",
                                 "-s", alignment,
                                 "-p", "1", "-f", "a", "-x", "1", "-#", str(num_reps),
                                 "-n", "BOOT{}".format(self.date)])
            else:
                subprocess.call([rax_ex, "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-s", alignment,
                                 "-p", "1", "-b", "1", "-#", str(num_reps),
                                 "-n", "BOOT{}".format(self.date)])
        outpath = "{}/RAxML_bootstrap.BOOT{}".format(self.rundir, self.date)
        return outpath



    def summarize_boot(self, besttreepath, bootpath, min_clade_freq=0.2):
        """Summarize the bootstrap proportions onto the ML tree"""
        summarized_tree_path = "{}/boot_summary.tre".format(self.rundir)
        sumtree_file = open(summarized_tree_path, 'w')
        subprocess.check_call(["sumtrees.py",
                               "-t", besttreepath,
                               "-f", str(min_clade_freq),
                               "-d0",
                               bootpath], stdout=sumtree_file)
        sumtree_file.close()
        debug(summarized_tree_path)
        return summarized_tree_path


    def calculate_final_tree(self, boot_reps=100):
        """Calculates the final tree using a trimmed alignment.

        :return: final PS data
        """
        debug("calculate final tree")
        debug("current alignment length {}".format(len(self.data.aln)))
        besttreepath = self.est_full_tree()
        if besttreepath is None:
            return None
        bootpath = self.calculate_bootstrap(num_reps=boot_reps)
        sumtreepath = self.summarize_boot(besttreepath, bootpath)
        self.replace_tre(sumtreepath, schema="nexus")
        self.data.write_files(direc=self.outputsdir)
        self.data.write_labelled(filename='updated_taxonname', label='^ot:ottTaxonName', direc=self.outputsdir)
        return 1


    def remove_blocklistitem(self):
        """This removes items from aln, and tree, if the corresponding Genbank identifer were added to the blocklist.

        Note, that seq that were not added because they were similar to the one being removed here, are lost
        (that should not be a major issue though, as in a new blast_run, new seqs from the taxon can be added.)
        """
        for tax in self.data.aln.taxon_namespace:
            gi_id = self.data.otu_dict[tax.label].get("^ncbi:gi")
            acc = self.data.otu_dict[tax.label].get("^ncbi:accession")
            if gi_id in self.blocklist or acc in self.blocklist:
                self.data.remove_taxa_aln_tre(tax.label)
                self.data.otu_dict[tax.label]['^physcraper:status'] = "deleted, Genbank identifier is part of blocklist"
        # this should not need to happen here: prune_short; instead...
        self.data.check_tre_in_aln()
        # self.data.prune_short()
        # debug(self.data.tre.as_string(schema='newick'))
