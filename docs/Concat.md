
`physcraper.Concat(object)`
=====================================================

Combines several physcraper runs into a concatenated alignment and calculates a phylogeny.
    
There are two options available, either data will be concatenated by random (per taxon name) or the 
    user provides a file which say, which sequences shall be concatenated.

User need to make sure, that there are at least some overlapping lineages.
    Do not concatenate data from the same loci (if you want to expand an alignment, run physcraper!).
    
####To build the class the following is needed:
  * **workdir_comb**: the path to your directory where the data shall be stored
  * **email**: your email address, currently used to retrieve missing taxon information

####During the initializing process the following self objects are generated:
  * **self.workdir**: the path to your directory
  * **self.sp_gi_comb**: dictonary - similar to otu_dcit of Physcraper class
      * key: taxon name
      * value: dictionary:
          * key: gene name
          * value: dictionary:
              * key: unique identifier of concat_class
              * value: dictionary - key-value-pairs:
                  * "gi_id": gi_id
                  * "seq": corresponding sequence
                  * "spn": taxon name
                  * "original_PS_id": otu_id from single-gene run
                  * "concat:status": "single run"/ "concatenated"
                  * "new tipname": taxon name plus number if there are more than a single
                                concatenated sequence per taxon
  * **self.single_runs**: dictionary
      * key: gene name, as provided by user in input
      * value: file containing the single gene Physcraper run, loaded from pickle
  * **self.sp_counter**: dictionary
      * key: taxon name
      * value: dictionary
          * key: gene name
          * value: number of sequences present for the given gene and taxon
  * **self.email**: email
  * **self.comb_seq**: dictionary
      * key: gene name
      * value: dictionary:
          * key: taxon name
          * value: sequence
  * **self.comb_gi**: dictonary # !!! Note can be easily combined with comb_seq
      * key: gene name
      * value: dictionary:
          * key: taxon name
          * value: unique identifier of concat_class
  * **self.aln_all**: dictionary
      * key: numbers from 1 to amount of loci concatenated
      * value: dendropy aln including empty seq
  * **self.num_of_genes**: number corresponding to the number of genes that shall be concatenated
  * **self.genes_present**: list of the gene names that shall be concatenated
  * **self.tre_as_start**: phylogeny used as starting tree for concatenated run, is the one with most tips present
  * **self.tre_start_gene**: corresponding name to tre_as_start
  * **self.short_concat_seq**: list of taxa that have only few sequences information/short genes
  * **self.concat_tips**: dictonary
      * key: otu.label from self.tre_as_start.taxon_namespace
      * value: taxon name
  * **self.concatfile**: path to file if user supplied concatenation file is used
  * **self.concatenated_aln**: concatenated alignment
  * **self.tmp_dict**: subset of self.sp_gi_comb
  * **self.part_len**: holds sequence partition position to write the partitioning file # ! TODO MK: might not need to be self
