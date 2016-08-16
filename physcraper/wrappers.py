#!/usr/bin/env python
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts,  PhyscraperScrape
from dendropy import DnaCharacterMatrix
import pickle
import sys
import os


def standard_run(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi):
    if os.path.isfile("{}/scrape.p".format(workdir)): 
        sys.stdout.write("Readloading from pickled scrapefile")
        scraper = pickle.load(open("{}/scrape.p".format(workdir),'rb'))
        scraper.repeat = 1
    else: 
            sys.stdout.write("setting up Data Object\n")
            sys.stdout.flush()
            #read the config file into a configuration object
            conf = ConfigObj(configfi)
            aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

            #Generate an linked Alignment-Tree-Taxa object
            data_obj = generate_ATT_from_phylesystem(aln=aln,
                                 workdir=workdir,
                                 study_id = study_id,
                                 tree_id = tree_id,
                                 phylesystem_loc = conf.phylesystem_loc)




            #Prune sequnces below a certain length threshold
            #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
            data_obj.prune_short()

            data_obj.write_files()
            data_obj.write_labelled()


            #Mapping identifiers between OpenTree and NCBI requires and identifier dict object
            ids = IdDicts(conf, workdir="example")


            #Now combine the data, the ids, and the configuration into a single physcraper scrape object
            scraper =  PhyscraperScrape(data_obj, ids, conf)
            #run the ananlyses
            scraper.run_blast()
            scraper.read_blast()
            scraper.remove_identical_seqs()
            scraper.generate_streamed_alignment()
    while scraper.repeat == 1: 
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()
