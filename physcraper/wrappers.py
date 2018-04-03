#!/usr/bin/env python
import sys
import os
import subprocess
import jsonpickle
import pickle
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from dendropy import DnaCharacterMatrix


def sync_ncbi(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["rsync", "av", "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz", "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp)])
    subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


def sync_ott(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["process_ott.sh", "".format(conf.ott_ncbi)])

def standard_run(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi):
    '''looks for a json file to continue run, or builds and runs
    new analysis for as long as new seqs are found'''
    conf = ConfigObj(configfi)
#    if os.path.isfile("{}/att_checkpoint.json".format(workdir)):
#        sys.stdout.write("Reloading data object from json scrapefile\n")
#        thawed = open("{}/att_checkpoint.json".format(workdir), 'r').readlines()
#        data_obj = jsonpickle.decode(thawed)
#        scraper.repeat = 1
    if os.path.isfile("{}/att_checkpoint.p".format(workdir)):
        sys.stdout.write("Reloading data object from pickle file\n")
        data_obj = pickle.load( open( "{}/att_checkpoint.p".format(workdir), "rb" ) )
#        scraper.repeat = 1
    else:
#            sync_names()
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
        #Mapping identifiers between OpenTree and NCBI requires and identifier dict object
        ids = IdDicts(conf, workdir="example")
        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label=  '^ot:ottTaxonName')
        data_obj.write_otus('otu_info', schema='table')
        data_obj.dump()
        #Mapping identifiers between OpenTree and NCBI requires and identifier dict object
    if os.path.isfile(conf.id_pickle):
        sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
#        thawed_id = open(conf.id_json, 'r').readlines()
#        ids = jsonpickle.decode(thawed_id)
#        scraper.repeat = 1
        ids = pickle.load(open(conf.id_pickle, "rb" ))
    else:
        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        ids.dump()
#Now combine the data, the ids, and the configuration into a single physcraper scrape object
    scraper = PhyscraperScrape(data_obj, ids, conf)
    #run the ananlyses
    scraper.run_blast()
    scraper.read_blast()
    scraper.remove_identical_seqs()
    scraper.generate_streamed_alignment()
    while scraper.repeat == 1:
        scraper.data.write_labelled()
        scraper.data.write_otus("otu_info", schema='table')
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()

