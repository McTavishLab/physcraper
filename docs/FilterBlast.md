`physcraper.FilterBlast(PhyscraperScrape)`
=====================================================

Takes the Physcraper Superclass and filters the ncbi blast results to only include a subset of the sequences.

They can be filtered by number or by rank and number. The feature can be useful for non-population-level studies,
e.g. analyses which require a single representative per taxon (threshold = 1) or to check the monophyly of taxa
without having to deal with over-representation of few taxa (e.g. threshold = 4, which allows to get a good overview
of what is available without having some taxa being represented by high numbers of sequences).
The second option (downtorank), which is optional, allows to filter according to taxonomic levels,
e.g. getting a number of representative sequences for a genus or lineage; it can also be used to deal with subspecies.


#### existing self objects are:
  * **self.sp_d**: dictionary
    * key = species name
    * value = dictionary:
        * key = otuID
        * value = otu_dict entry
  * **self.**: dictionary
    * key = species name
    * value = dictionary (Is overwritten every 'round')
        * key = otuID
        * value = seq.
  * **self.filtered_seq**: dictionary. Is used as the self.new_seqs equivalent from Physcraper, just with fewer seqs.
                    Is overwritten every 'round'
    * key = otuID,
    * val = seq.
  * **self.downtorank**: optional string defining the level of taxonomic filtering, e.g. "species", "genus"
  