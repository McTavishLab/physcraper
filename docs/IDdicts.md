
`physcraper.IdDicts(object)`
=====================================================


Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

####To build the class the following is needed:
    
  * **config_obj**: Object of class config (see above)
  * **workdir**: the path to the assigned working directory
  * **mrca**: mrca as defined by input, can be a class

####During the initializing process the following self objects are generated:
  * **self.workdir**: contains path of working directory
  * **self.config**: contains the Config class object
  * **self.ott_to_ncbi**: dictionary
  
      * key: OToL taxon identifier
      * value: ncbi taxon identifier
  * **self.ncbi_to_ott**: dictionary
  
      * key: OToL taxon identifier
      * value: ncbi taxon identifier
  * **self.ott_to_name**: dictionary
  
      * key: OToL taxon identifier
      * value: OToL taxon name
  * **self.acc_ncbi_dict**: dictionary
  
      * key: Genbank identifier
      * value: ncbi taxon identifier
  * **self.spn_to_ncbiid**: dictionary
  
      * key: OToL taxon name
      * value: ncbi taxon identifier
  * **self.ncbiid_to_spn**: dictionary
  
      * key: ncbi taxon identifier
      * value: ncbi taxon name

  * **self.mrca_ott**: user defined list of mrca OTT-ID's
  * **self.mrca_ncbi**: set, which is fed by self.get_ncbi_mrca()

  * **Optional**:
      * depending on blasting method:
       * self.ncbi_parser: for local blast, initializes the ncbi_parser class, that contains information about rank and identifiers
       * self.otu_rank: for remote blast to store the rank information
