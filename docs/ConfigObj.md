`physcraper.ConfigObj(object)`
=====================================================
Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store.

####To build the class the following is needed:
  * **configfi**: a configuration file in a specific format,
   
        e.g. to read in self.e_value_thresh.
            The file needs to have a heading of the format: [blast] and then somewhere below that heading
            a string e_value_thresh = value

  * **interactive**: defaults to True, is used to interactively update the local blast databases

####During the initializing process the following self objects are generated:
  * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out:
                            https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
  * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
  * **self.minlen**: value from 0 to 1. Defines how much shorter new seq can be compared to input
  * **self.get_ncbi_taxonomy**: Path to sh file doing something...
  * **self.ncbi_dmp**: path to file that has gi numbers and the corresponding ncbi tax id's
  * **self.phylesystem_loc**: defines which phylesystem for OpenTree datastore is used. 
      The default is api, but can run on local version too. 
  * **self.ott_ncbi**: file containing OTT id, ncbi and taxon name (??)
  * **self.id_pickle**: path to pickle file
  * **self.email**: email address used for blast queries
  * **self.blast_loc**: defines which blasting method to use:
      * either web-query (=remote)
      * from a local blast database (=local)
  * **self.num_threads**: number of cores to be used during a run
  * **self.url_base**: 
      * if blastloc == remote: it defines the url for the blast queries.
      * if blastloc == local: url_base = None
  * **self.unmapped**: used for OToL original tips that can not be assigned to a taxon
      * keep: keep the unmapped taxa and asign them to life
      * remove: remove the unmapped taxa from aln and tre
  * **optional self.objects**:
      * if blastloc == local:
          * self.blastdb: this defines the path to the local blast database
          * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
          * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
     