
`physcraper.AlignTreeTax(object)`
=====================================================

wrap up the key parts together, requires OTT_id, and names must already match.
Hypothetically, all the keys in the  otu_dict should be clean.

####To build the class the following is needed:

  * **newick**: dendropy.tre.as_string(schema=schema_trf) object
  * **otu_dict**: json file including the otu_dict information generated earlier
  * **alignment**: dendropy DNACharacterMatrix object
  * **ingroup_mrca**: OToL identifier of the group of interest, either subclade as defined by user or of
                    all tiplabels in the phylogeny
  * **workdir**: the path to the corresponding working directory
  * **schema**: optional argument to define tre file schema, if different from "newick"

####During the initializing process the following self objects are generated:

  * **self.aln**: contains the alignment and which will be updated during the run
  * **self.tre**: contains the phylogeny, which will be updated during the run
  * **self.otu_dict**: dictionary with taxon information and physcraper relevant stuff
        
       * key: a unique identifier (otu plus either "tiplabel of phylogeny" or for newly found sequences
            PS_number.
       * value: dictionary with the following key:values:
              
            * '^ncbi:gi': GenBank identifier - deprecated by Genbank - only older sequences will have it
            * '^ncbi:accession': Genbanks accession number
            * '^ncbi:title': title of Genbank sequence submission
            * '^ncbi:taxon': ncbi taxon identifier
            * '^ot:ottId': OToL taxon identifier
            * '^physcraper:status': contains information if it was 'original', 'queried', 'removed',
                                'added during filtering process'
            * '^ot:ottTaxonName': OToL taxon name
            * '^physcraper:last_blasted': contains the date when the sequence was blasted.
                                        
                 If the year is different from the 20th century, it tells us
                 something about the initial status:
                 * 1800 = never blasted, not yet considered to be added
                 * 1900 = never blasted and not added - see status for more information
                 * this century = blasted and added.
            * '^user:TaxonName': optional, user given label from OtuJsonDict
            * "^ot:originalLabel" optional, user given tip label of phylogeny
  * **self.ps_otu**: iterator for new otu IDs, is used as key for self.otu_dict
  * **self.workdir**: contains the path to the working directory, if folder does not exists it is generated.
  * **self.ott_mrca**: OToL taxon Id for the most recent common ancestor of the ingroup
  * **self.orig_seqlen**: list of the original sequence length of the input data
  * **self.gi_dict**: dictionary, that has all information from sequences found during the blasting.
        
    * key: GenBank sequence identifier
    * value: dictionary, content depends on blast option, differs between webquery and local blast queries
        * keys - value pairs for local blast:
            * '^ncbi:gi': GenBank sequence identifier
            * 'accession': GenBank accession number
            * 'staxids': Taxon identifier
            * 'sscinames': Taxon species name
            * 'pident': Blast  percentage of identical matches
            * 'evalue': Blast e-value
            * 'bitscore': Blast bitscore, used for FilterBlast
            * 'sseq': corresponding sequence
            * 'title': title of Genbank sequence submission
        * key - values for web-query:
            * 'accession':Genbank accession number
            * 'length': length of sequence
            * 'title': string combination of hit_id and hit_def
            * 'hit_id': string combination of gi id and accession number
            * 'hsps': Bio.Blast.Record.HSP object
            * 'hit_def': title from GenBank sequence
        * optional key - value pairs for unpublished option:
            * 'localID': local sequence identifier
  * **self._reconciled**: True/False,
  * **self.unpubl_otu_json**: optional, will contain the OTU-dict for unpublished data, if that option is used

####Following functions are called during the init-process:

  * **self._reconcile_names()**: 
        removes taxa, that are not found in both, the phylogeny and the aln and changes their names????

####The physcraper class is then updating: 
  * self.aln, self.tre and self.otu_dict, self.ps_otu, self.gi_dict


