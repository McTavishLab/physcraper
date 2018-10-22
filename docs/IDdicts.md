
`physcraper.IdDicts(object)`
=====================================================


Wraps up the annoying conversions. Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

####To build the class the following is needed:
    
  * **config_obj**: Object of class config (see above)
  * **workdir**: the path to the assigned working directory

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
  * **self.gi_ncbi_dict**: dictionary
  
      * key: Genbank identifier
      * value: ncbi taxon identifier
  * **self.spn_to_ncbiid**: dictionary
  
      * key: OToL taxon name
      * value: ncbi taxon identifier
  * **self.ncbiid_to_spn**: dictionary
  
      * key: ncbi taxon identifier
      * value: ncbi taxon name
  * **Optional**:
        
       * self.ncbi_parser: initializes the ncbi_parser class, that contains information about rank and identifiers
