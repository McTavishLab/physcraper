import random
from copy import deepcopy
from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel


from physcraper import PhyscraperScrape, filter_by_local_blast
from physcraper import ncbi_data_parser

from physcraper.helpers import cd


_DEBUG = 1
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)

class FilterBlast(PhyscraperScrape):
    """Takes the Physcraper Superclass and filters the ncbi blast results to only include a subset of the sequences.

    They can be filtered by number or by rank and number. The feature can be useful for non-population-level studies,
    e.g. analyses which require a single representative per taxon (threshold = 1) or to check the monophyly of taxa
    without having to deal with over-representation of few taxa (e.g. threshold = 4, which allows to get a good overview
    of what is available without having some taxa being represented by high numbers of sequences).
    The second option (downtorank), which is optional, allows to filter according to taxonomic levels, e.g. getting
    a number of representative sequences for a genus or lineage it can also be used to deal with subspecies.

    Existing self objects are:

        self.sp_d: dictionary

                key = species name/id
                value = list with otu_dict entry

        self.sp_seq_d: dictionary

                key = species name/id
                value = dictionary (Is overwritten every 'round')

                    key = otuID
                    value = seq.
        self.filtered_seq: dictionary. Is used as the self.new_seqs equivalent from Physcraper, just with fewer seqs. Is overwritten every 'round'

                key = otuID,
                val = seq.
        self.downtorank: optional string defining the level of taxonomic filtering, e.g. "species", "genus"
    """

    # TODO MK: self.sp_d = {} does not need to contain all otu_dict info, key is sufficient

    def __init__(self, data_obj, ids, ingroup_mrca=None):
        super(FilterBlast, self).__init__(data_obj, ids, ingroup_mrca)
        debug("start derived class init")
        # additional things that are needed for the filtering process
        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}
        self.downtorank = None
        self.threshold = None

    def add_setting_to_self(self, downtorank, threshold):
        """
        Add FilterBlast items to self.

        Currently used by some wrapper functions.

        :param downtorank: rank which defines your level of Filtering
        :param threshold: number, defining how many seq per rank do you want to keep
        :return:
        """
        self.threshold = threshold
        self.downtorank = downtorank

    def sp_dict(self, downtorank=None):
        """Takes the information from the Physcraper otu_dict and makes a dict with species id as key and
        the corresponding seq information from aln and blast seq, it returns self.sp_d.

        This is generated to make information for the filtering class more easily available. self.sp_d sums up which
        information are available per taxonomic concept and which have not already been removed during either
        the remove_identical_seq steps or during a filtering run of an earlier cycle.

        Note: has test, test_sp_d.py

        :param downtorank: string defining the level of taxonomic filtering, e.g. "species", "genus"
        :return: self.sp_d
        """
        #Edited to pass in list of otu_ids rather than set of dictionaries, to make getting squence by id easier in sp_seq_d
        self.downtorank = downtorank
        debug("make sp_dict")
        self.sp_d = {}
        for otu_id in self.data.otu_dict:
            if self.data.otu_dict[otu_id]['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                tax_id = self.data.otu_dict[otu_id].get('^ncbi:taxon') 
                assert tax_id not in set([0, None]) # every OTU must have a taxon_id for filter blast
                    # we cannot include unmapped taxa in fliter blast.
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
                if tax_id in self.sp_d:
                    self.sp_d[tax_id].append(otu_id)
                else:
                    self.sp_d[tax_id] = [otu_id]
        return self.sp_d

    def make_sp_seq_dict(self):
        """Uses the sp_d to make a dict with species names as key1, key2 is gb_id/sp.name and value is seq

        This is used to select representatives during the filtering step, where it selects how many sequences per
        species to keep in the alignment. It will only contain sp that were not removed in an earlier cycle of the
        program.

        Note: has test, test_sp_seq_d.py

        return: self.sp_seq_d
        """
        debug("make_sp_seq_dict")
        for key in self.sp_d:
            # loop to populate dict. key1 = sp id, key2= gb id, value = seq,
            # number of items in key2 will be filtered according to threshold and already present seq
            seq_d = {}
            for otu_id in self.sp_d[taxon]:
#                debug("otu id is {}".format(otu_id))
                try:
                    seq = self.data.aln[otu_id].symbols_as_string()
                except KeyError:
                    assert otu in self.new_seqs_otu_id # if it is not already in the alignment it must be in new seqs
                    seq = self.new_seqs_otu_id[otu]
                seq = seq.replace("-", "")
                seq = seq.replace("?", "")
                seq_d[otu_id] = seq
            self.sp_seq_d[taxon] = seq_d
        return

#    def select_seq_by_local_blast(self, seq_d, fn, count):
        # """
        # Selects number of sequences from local_blast to fill up sequences to the threshold.
        # Count is the return value from self.count_num_seq(tax_id)["seq_present"], that tells the program
        # how many sequences for the taxon are already available in the aln.

        # It will only include species which have a blast score of mean plus/minus sd.
        # Uses the information returned by read_local_blast_query()
        # to select which sequences shall be added in a filtered run.

        # Note: has test,test_select_seq_by_local_blast.py

        # :param seq_d: is the value of self.sp_d (= another dict)
        # :param fn: refers to a filename to find the local blast file produced before,
        #             which needs to be read in by read_local_blast_query() - currently tax_id
        # :param count: self.count_num_seq(tax_id)["seq_present"]
        # :return: self.filtered_seq
        # """
        # debug("select_seq_by_local_blast")
        # no_similar_seqs = 0
        # seq_blast_score = filter_by_local_blast.read_filter_blast(self.workdir, seq_d, fn)
        # random_seq_ofsp = {}

        # if seq_blast_score != {}:
        #     if (self.threshold - count) <= 0:
        #         debug("already too many samples of sp in aln, skip adding more.")
        #     # exact amount of seq present which need to be added
        #     elif len(seq_blast_score.keys()) == (self.threshold - count):
        #         random_seq_ofsp = seq_blast_score
        #     elif len(seq_blast_score.keys()) > (
        #             self.threshold - count):  # more seq available than need to be added, choose by random
        #         debug("choose random")
        #         random_seq_ofsp = random.sample(seq_blast_score.items(), (self.threshold - count))
        #         random_seq_ofsp = dict(random_seq_ofsp)
        #     elif len(seq_blast_score.keys()) < (
        #             self.threshold - count):  # less seq available than need to be added, just use all
        #         debug("add all")
        #         # print(len(seq_blast_score.keys()))
        #         # print(self.threshold - count)
        #         random_seq_ofsp = seq_blast_score
        # # no similar seq found. think about what to do. was the other seq already present?
        # else:
        #     debug("blast did not find similar seqs")
        #     if len(seq_d.keys()) > 2 and no_similar_seqs == 0:  # try with different seq to blast
        #         debug("blast with different seq...")
        #         # all the next line is from how_many_seq_to_keep()
        #         blast_seq_id = seq_d.keys()[1]  # seq 1 instead of 0 now
        #         seq = seq_d[blast_seq_id]
        #         filter_by_local_blast.write_filterblast_query(self.workdir, blast_seq_id, seq,
        #                                                       fn=fn)  # blast guy
        #         blast_db = seq_d.keys()[2:]
        #         for blast_key in blast_db:
        #             seq = seq_d[blast_key]
        #             filter_by_local_blast.write_filterblast_db(self.workdir, blast_key, seq,
        #                                                           fn=fn)
        #         # make local blast of sequences
        #         filter_by_local_blast.run_filter_blast(self.workdir, fn, fn)
        #         no_similar_seqs = 1
        #         count_dict = self.count_num_seq(fn)
        #         seq_present = count_dict["seq_present"]
        #         if len(seq_d) + seq_present >= self.threshold:
        #             self.select_seq_by_local_blast(seq_d, fn, count)
        #         elif len(seq_d) + seq_present < self.threshold:
        #             self.add_all(fn)
        #     elif len(seq_d.keys()) > 2 and no_similar_seqs == 1:  # also with different seq no result, select random seq!
        #         if len(seq_d.keys()) == (self.threshold - count):  # exact amount of seq present which need to be added
        #             for item in seq_d.keys():
        #                 random_seq_ofsp[item] = seq_d[item]
        #         elif len(seq_d.keys()) > (
        #                 self.threshold - count):  # more seq available than need to be added, choose by random
        #             debug("choose random - else")
        #             random_seq = random.sample(seq_d.items(), (self.threshold - count))
        #             for item in random_seq.keys():
        #                 random_seq_ofsp[item] = random_seq[item]
        #         elif len(seq_blast_score.keys()) < (
        #                 self.threshold - count):  # less seq available than need to be added, just use all
        #             debug("add all - else")
        #             for item in seq_d.keys():
        #                 random_seq_ofsp[item] = seq_d[item]
        #     else:
        #         seq_id = seq_d.keys()[1]
        #         seq = seq_d[seq_id]
        #         random_seq_ofsp[seq_id] = seq
        # # debug(random_seq_ofsp)
        # if len(random_seq_ofsp) > 0:  # add everything to filtered seq
        #     for key, val in random_seq_ofsp.items():
        #         self.filtered_seq[key] = val
        # return self.filtered_seq

    def select_seq_by_length(self, taxon_id, count):
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
        max_len = max(self.sp_seq_d[taxon_id].values())

        seq_w_maxlen = {}
        for key, val in self.sp_seq_d[taxon_id].items():
            for otu_id in self.sp_d[taxon_id]:
                item = self.data.otu_dict[otu_id]
                if '^ncbi:accession' in item and item['^ncbi:accession'] == key:
                    if item['^physcraper:status'].split(' ')[0] != ["added", "deleted", "original", "new"]:
                        if len(val) == len(max_len):
                                seq_w_maxlen[key] = val
        if (self.threshold - count) <= 0:
            debug("already to many samples of sp in aln, skip adding more.")
            random_seq_ofsp = None
        elif len(seq_w_maxlen) == (self.threshold - count):
            random_seq_ofsp = seq_w_maxlen
        elif len(seq_w_maxlen) > (self.threshold - count):
            random_seq_ofsp = random.sample(seq_w_maxlen.items(), (self.threshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        else:
            toselect = range(len(seq_w_maxlen), (self.threshold - count))
            keymax = seq_w_maxlen.keys()
            subdict = {k: v for k, v in self.sp_seq_d[taxon_id].items() if k not in keymax}
            second_len = max(subdict.values())
            seq2len = {}
            for key, val in subdict.items():
                if len(val) == len(second_len):
                    seq2len[key] = val
            random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
            random_seq_ofsp = dict(random_seq_ofsp)
            random_seq_ofsp.update(seq_w_maxlen)
        if (self.threshold - count) >= 1:
            assert random_seq_ofsp is not None
        if random_seq_ofsp is not None:
            for key in random_seq_ofsp.keys():
                self.filtered_seq[key] = random_seq_ofsp[key]

    def add_all(self, key):
        """It adds all seq to filtered_seq dict as the number of sequences present is smaller than the threshold value.

        It is only used, when sequences selection happens via blasting.

        Note: has test, test_add_all.py

        :param key: key of self.sp_d (taxon id)
        :return: self.filtered_seq
        """
        debug('add_all')
        for otu_id in self.sp_d[key]:
            otu_info = self.data.otu_dict[otu_id]
            if '^physcraper:status' in otu_info:
                # debug(otu_id)
                if otu_info['^physcraper:status'].split(' ')[0] not in self.seq_filter:
                    if otu_info['^physcraper:last_blasted'] == None \
                            and otu_info['^physcraper:status'] != "original":
                        seq = self.new_seqs_otu_id[otu_id]
                        self.filtered_seq[otu_id] = seq
        return self.filtered_seq

#    def loop_for_write_blast_files(self, tax_id):
        # """This loop is needed to be able to write the local blast files for the filtering step correctly.

        # Function returns a filename for the filter blast.

        # Note: has test,test_loop_for_blast.py

        # :param key: key of self.sp_d (taxon id)
        # :return: name of the blast file
        # """
        # debug("loop_for_write_blast_files")
        # aln_otus = set([taxon.label for taxon in self.data.aln])
        # query_otu = None
        # db_otus = []
        # for otu_id in self.sp_d[tax_id]: 
        #     otu_info = self.data.otu_dict[otu_id]
        #     if otu_id in aln_otus:
        #         query_otu = otu_id #we end up overwriting the query file repeatedly. Might as well just choose one otu and write it once.
        #     elif '^physcraper:status' in otu_info and otu_info['^physcraper:status'].split(' ')[0] not in self.seq_filter: # these are the new sequences that haven't been filtered out
        #         db_otus.append(otu_id)
        #         assert otu_id in self.new_seqs_otu_id.keys()
        #         seq = self.new_seqs_otu_id[otu_id]
        #         filter_by_local_blast.write_filterblast_db(self.workdir, gb_id, seq, fn=tax_id)
        #     else:
        #         debug("otu_id {} was not in the alignemnt, but was filtered out due to {}".format(otu_id, otu_info['^physcraper:status']))
        # assert query_otu is not None#at least 1 otu must be in the alignment
        # debug("for taxon {} will use otu {} for query".format(tax_id, query_otu))
        # query_seq = self.data.aln[query_otu]
        # filter_by_local_blast.write_filterblast_query(self.workdir, query_otu, query_seq, fn=tax_id)
        # return tax_id


    def count_num_seq(self, tax_id):
        """Counts how many sequences there are for a tax_name, excluding sequences that have not been added
        during earlier cycles.

        Function is only used in how_many_sp_to_keep().

        :param tax_id: key from self.sp_seq_d
        :return: dict which contains information of how many seq are already present in aln, how many new seq have been
                found and if the taxon is a new taxon or if seq are already present
        """
        debug("count_num_seq for tax_id {}".format(tax_id))
        seq_added = 0
        original = 0
        new_taxon = True
        query_count = 0
        seq_in_aln = 0
        for otu_id in self.sp_d[tax_id]:
            item = self.data.otu_dict[otu_id]
            aln_otus = set([taxon.label for taxon in self.data.aln])
            if otu_id in aln_otus:
                seq_in_aln += 1
                new_taxon = False
            status = item.get('^physcraper:status')
            assert status is not None
            if status.split(' ')[0] not in self.seq_filter:
#                debug(item['^physcraper:status'])
                item_split = item['^physcraper:status'].split(' ')[0]
                if item["^physcraper:status"] == "query" or item_split == "new" or item_split == "added,":
                    query_count += 1
                if item["^physcraper:status"] == 'added as representative of taxon':
                    seq_added += 1
                    new_taxon = False
                if item_split == "original":
                    original += 1
                    new_taxon = False
        seq_present = seq_added + original
        assert seq_in_aln == seq_present
        # if item_split == "added," or item_split == "original":
        count_dict = {
            "seq_present": seq_added + original,
            "query_count": query_count,
            "new_taxon": new_taxon,
        }
        if new_taxon is False:
            assert original != 0 or seq_added != 0, ("count_dict `%s` has more seq added than threshold: 0." % count_dict)
        if new_taxon is True:
            assert original == 0, ("count_dict `%s` has more original seq than allowed for new taxon." % count_dict)
            assert seq_added == 0, ("count_dict `%s` has more seq added than allowed for new taxon." % count_dict)
        # debug([seq_added, original, self.threshold])
        if original < self.threshold:
            assert seq_added <= self.threshold, ("count_dict `%s` has more seq added than threshold." % count_dict)
        elif original > self.threshold:
            sys.stdout.write("already more originals than requested by threshold...\n")
        else:
            assert seq_added + original <= self.threshold, \
                "seq_added ({}) and original ({}) have more than threshold ({}).".format(seq_added, original, self.threshold)
        return count_dict

    def how_many_sp_to_keep(self, selectby):
        """Uses the sp_seq_d and places the number of sequences according to threshold into the self.filterdseq_dict.

        This is essentially the key function of the Filter-class, it wraps up everything.

        :param selectby: mode of sequence selection, defined in input
        :return: nothing specific, it is the function, that completes the self.filtered_seq, which contains the filtered
                sequences that shall be added.
        """
        debug("how_many_sp_to_keep")
        # self.threshold = threshold
        for tax_id in self.sp_d:
            count_dict = self.count_num_seq(tax_id)
            seq_present = count_dict["seq_present"]
            query_count = count_dict["query_count"]
            new_taxon = count_dict["new_taxon"]
            #debug(count_dict)
            # debug(tax_id)
            if seq_present <= self.threshold:  # add seq to aln
                if seq_present + query_count <= self.threshold:  # to add all stuff to self.filtered_seq[gi_n]
                    debug("taxon {}, has {} seqs present, and we found {}, which totals less than the threshold {}, so we will add all".format(tax_id, seq_present, query_count, self.threshold))
                    self.add_all(tax_id)
                else:  # filter number of sequences
                    debug("taxon {}, has {} seqs present, and we found {}, which totals MORE than the threshold {}, so we will select by {}".format(tax_id, seq_present, query_count, self.threshold,selectby))
                    if tax_id in self.sp_seq_d.keys():
                        if selectby == "length":
                            self.select_seq_by_length(tax_id, seq_present)
                        elif selectby == "blast":
                            if seq_present == 0 and new_taxon is True and query_count >= 1:  # if new taxon
                                # debug("new taxon")
                                debug(self.sp_seq_d[tax_id].keys())
                                blast_seq_id = self.sp_seq_d[tax_id].keys()[0]
                                seq = self.sp_seq_d[tax_id][blast_seq_id]
                                filter_by_local_blast.write_filterblast_query(self.workdir, blast_seq_id, seq,
                                                                              fn=tax_id)  # blast guy
                                blast_db = self.sp_seq_d[tax_id].keys()[1:]
                                for blast_key in blast_db:
                                    seq = self.sp_seq_d[tax_id][blast_key]
                                    filter_by_local_blast.write_filterblast_db(self.workdir, blast_key, seq, fn=tax_id)
                                # make local blast of sequences
                                filter_by_local_blast.run_filter_blast(self.workdir, tax_id, tax_id)
                                debug(self.sp_seq_d[tax_id])
                                if len(self.sp_seq_d[tax_id]) + seq_present >= self.threshold:
                                    self.select_seq_by_local_blast(self.sp_seq_d[tax_id], tax_id, seq_present)
                                elif len(self.sp_seq_d[tax_id]) + seq_present < self.threshold:
                                    self.add_all(tax_id)
                            elif 1 <= seq_present < self.threshold and new_taxon is False and query_count != 0:
                                if query_count + seq_present > self.threshold:
                                    taxonfn = self.loop_for_write_blast_files(tax_id)
                                    if self.downtorank is not None:
                                        taxonfn = tax_id
                                    filter_by_local_blast.run_filter_blast(self.workdir, taxonfn, taxonfn)
                                    self.select_seq_by_local_blast(self.sp_seq_d[tax_id], taxonfn, seq_present)
                                elif query_count + seq_present <= self.threshold:
                                    self.add_all(tax_id)
        return

    def replace_new_seq(self):
        """Function to replace self.new_seqs and self.new_seqs_otu_id with the subset of filtered sequences.

        This is the final step in the FilterBlast class, from here it goes back to PhyScraper.

        :return: subsets of self.new_seqs and self.new_seqs_otu_id
        """
        debug("replace new seq")
        keylist = self.filtered_seq.keys()
        if not self.unpublished:
            keylist = [x for x in keylist if x[:6] != "unpubl"]
        seq_not_added = self.new_seqs_otu_id.keys()
        reduced_new_seqs_dic = {}
        for otu_id in seq_not_added:
            self.data.otu_dict[otu_id]['^physcraper:status'] = 'not added, there are enough seq per sp in tre'
        for otu_id in keylist:
                reduced_new_seqs_dic[otu_id] = self.filtered_seq[otu_id]
                self.data.otu_dict[otu_id]['^physcraper:status'] = 'added as representative of taxon'

        reduced_new_seqs = {k: self.filtered_seq[k] for k in keylist}
        with open(self.logfile, "a") as log:
            log.write("{} sequences added after filtering, of {} before filtering\n".format(len(reduced_new_seqs_dic),
                                                                                            len(self.new_seqs_otu_id)))
        self.new_seqs = deepcopy(reduced_new_seqs)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic)
        # set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        return


####################
#Funcs below here should be moved to a separate script at some point
