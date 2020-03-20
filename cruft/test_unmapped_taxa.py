import physcraper
import sys
import pickle
import os
from dendropy import DnaCharacterMatrix

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "ot_350"
tree_id = "Tr53296"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "tests/output/opentree_unmappedtaxa"
absworkdir = os.path.abspath(workdir)

configfi = "tests/data/test.config"

def test_0():
    
    sys.stdout.write("\n\n No files present\n\n")
    conf = physcraper.ConfigObj(configfi)
    conf.unmapped = 'keep'
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                workdir=workdir,
                                config_obj=conf,
                                study_id = study_id,
                                tree_id = tree_id,
                                phylesystem_loc = conf.phylesystem_loc)
    # physcraper.debug(len(data_obj.aln.taxon_namespace))
    data_obj.workdir = absworkdir

    ids =  physcraper.IdDicts(conf, workdir=workdir)
    # physcraper.debug(os.getcwd())
    data_obj.write_files()
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    # physcraper.debug(len(scraper.data.aln.taxon_namespace))
    # physcraper.debug("scraper obj made")
    num_keep = len(scraper.data.aln.taxon_namespace)
    # physcraper.debug(num_keep)




    conf = physcraper.ConfigObj(configfi, interactive=False)
    conf.unmapped = 'remove'


    data_obj.workdir = absworkdir

    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)


    scraper2 = physcraper.PhyscraperScrape(data_obj, ids)
    num_remove = len(scraper2.data.aln.taxon_namespace)
    dict_id = 0
    for tax in scraper.data.aln.taxon_namespace:
        if '^ot:ottId' in scraper.data.otu_dict[tax.label]:
            dict_id = dict_id + 1
    # print(num_remove, num_keep, dict_id)
    assert num_remove <= num_keep - 1
    assert num_keep == dict_id
