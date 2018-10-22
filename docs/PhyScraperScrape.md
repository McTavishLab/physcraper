`physcraper.PhyscraperScrape(object)`
=====================================================

This is the class that does the perpetual updating

####To build the class the following is needed:

  * **data_obj**: Object of class ATT (see above)
  * **ids_obj**: Object of class IdDict (see above)

####During the initializing process the following self.objects are generated:

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
  * **self.mrca_ncbi**: ncbi identifier of mrca

  * **self.tmpfi**: path to a file or folder???
  * **self.blast_subdir**: path to folder that contains the files writen during blast

  * **self.newseqs_file**: filename of files that contains the sequences from self.new_seqs_otu_id
  * **self.date**: Date of the run - may lag behind real date!
  * **self.repeat**: either 1 or 0, it is used to determine if we continue updating the tree, no new seqs found = 0
  * **self.newseqsgi**: list of all gi_ids that were passed into remove_identical_seq(). Used to speed up adding process
  * **self.blacklist**: list of gi_id of sequences that shall not be added or need to be removed. Supplied by user.
  * **self.gi_list_mrca**: list of all gi_ids available on GenBank for a given mrca. Used to limit possible seq to add.
  * **self.seq_filter**: list of words that may occur in otu_dict.status and which shall not be used
                in the building of FilterBlast.sp_d (that's the main function), but it is also used as assert
                statement to make sure unwanted seqs are not added.
  * **self.unpublished**: True/False. Used to look for local unpublished seq that shall be added if True.
  * **self.path_to_local_seq:** Usually False, contains path to unpublished sequences if option is used.

####Following functions are called during the init-process:
* **self.reset_markers()**: 
     adds things to self: I think they are used to make sure certain function run, if program crashed and pickle file is read in.
    * self._blasted: 0/1, if run_blast() was called, it is set to 1 for the round.
    * self._blast_read: 0/1, if read_blast() was called, it is set to 1 for the round.
    * self._identical_removed: 0
    * self._query_seqs_written: 0/1, if write_query_seqs() was called, it is set to 1 for the round.
    * self._query_seqs_aligned: 0
    * self._query_seqs_placed: 0/1, if place_query_seqs() was called, it is set to 1 for the round.
    * self._reconciled: 0
    * self._full_tree_est: 0/1, if est_full_tree() was called, it is set to 1 for the round.
  