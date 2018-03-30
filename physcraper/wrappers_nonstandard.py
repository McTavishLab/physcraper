#!/usr/bin/env python
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts,  PhyscraperScrape
from dendropy import DnaCharacterMatrix
import pickle
import sys
import os
import subprocess
import csv
from ete2 import NCBITaxa


def sync_ncbi(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["rsync", "av", "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz", "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp)])
    subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


def sync_ott(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["process_ott.sh",  "".format(conf.ott_ncbi)])


##generates IdDicts physcrapper class        
def get_ottid(configfi, cwd):
                     conf = ConfigObj(configfi)
                     ids = IdDicts(conf, cwd)  
                     return(ids)            


print("load OtuJsonDict func")    

def OtuJsonDict(id_to_spn, configfi):
    """Make otu json dict, which is also produces within the openTreeLife-query"""
    cwd = os.getcwd()  
    ## reads input file into the var spInfo
    with open(id_to_spn, mode='r') as idtospn:
        reader = csv.reader(idtospn)
        spInfo = dict((rows[0], rows[1]) for rows in reader)
    #print(spInfoDict) 
 
    ###generate spinfodict
    
    ottdic = get_ottid(configfi, cwd) 
    print(ottdic)
    ncbi = NCBITaxa()    
    
    spInfoDict = {}
    for item in spInfo:
        spn = spInfo[item].replace("_", " ")
        name2taxid = ncbi.get_name_translator([spn])
       
        
        if len(name2taxid.items())>=1:
            ncbiid = name2taxid.items()[0][1][0]
            ott = ottdic.ncbi_to_ott[ncbiid]
            spn = ottdic.ott_to_name[ott]
            get_info = {'^ncbiID': ncbiid, '^ot:ottTaxonName': spn, 'ot:ottId': ott, '^user:TaxonName': spInfo[item],  '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01" }  
            spInfoDict[item] = get_info
        else:
            
            spInfoDict[item] = {'^user:TaxonName': spInfo[item],  '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01"}
    print(spInfoDict)
    return  spInfoDict 

    

def own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 spInfoDict,
                 configfi):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found'''
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile\n")
        scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        scraper.repeat = 1
    elif os.path.isfile("{}/ATT_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile\n")
        scraper = pickle.load(open("{}/ATT_checkpoint.p".format(workdir),'rb'))
        scraper.repeat = 1    
    else:   
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        print(seqaln, mattype)
        #aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
        
        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             treefile=trfn,
                             schema_trf = schema_trf,
                             otu_json=spInfoDict,
                             ingroup_mrca=None)

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()

        data_obj.write_files()

        data_obj.write_labelled( label='user:TaxonName')
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()
        
        ids = IdDicts(conf, workdir="example")


        #Now combine the data, the ids, and the configuration into a single physcraper scrape object
        scraper =  PhyscraperScrape(data_obj, ids, conf)
        #run the ananlyses
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
        scraper.dump()
        scraper.how_many_sp_to_keep(treshold=treshhold)

        scraper.generate_streamed_alignment()
    while scraper.repeat == 1: 
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
#        scraper.how_many_sp_to_keep(treshold=treshhold)

        scraper.generate_streamed_alignment()


# ott_ids = get_subtree_otus(nexson,
#                                tree_id=tree_id,
#                                subtree_id="ingroup",
#                                return_format="ottid")
#     ott_mrca = get_mrca_ott(ott_ids)
#     newick = extract_tree(nexson,
#                           tree_id,
#                           PhyloSchema('newick',
#                                       output_nexml2json='1.2.1',
#                                       content="tree",
#                                       tip_label="ot:originalLabel"))
#     newick = newick.replace(" ", "_") #UGH Very heavy handed, need to make sure happens on alignement side as well.
#     tre = Tree.get(data=newick,
#                    schema="newick",
#                    preserve_underscores=True,
#                    taxon_namespace=aln.taxon_namespace)
#     otus = get_subtree_otus(nexson, tree_id=tree_id)
#     otu_dict = {}
#     orig_lab_to_otu = {}
#     treed_taxa = {}
#     for otu_id in otus:
#         otu_dict[otu_id] = extract_otu_nexson(nexson, otu_id)[otu_id]
#         otu_dict[otu_id]['^physcraper:status'] = "original"
#         otu_dict[otu_id]['^physcraper:last_blasted'] = "1900/01/01"
#         orig = otu_dict[otu_id].get(u'^ot:originalLabel').replace(" ", "_")
#         orig_lab_to_otu[orig] = otu_id
#         treed_taxa[orig] = otu_dict[otu_id].get(u'^ot:ottId')
#     for tax in aln.taxon_namespace:
#         try:
#             tax.label = orig_lab_to_otu[tax.label].encode('ascii')
#         except KeyError:
#             sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignement. This may indicate a mismatch between tree and alignement\n".format(tax.label))
#    #need to prune tree to seqs and seqs to tree...     
