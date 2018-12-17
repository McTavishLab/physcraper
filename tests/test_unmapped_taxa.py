import physcraper
import sys
import pickle
import os
from dendropy import DnaCharacterMatrix

#Use OpenTree phylesystem identifiers to get study and tree
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype = "fasta"
workdir = "tests/output/opentree_unmappedtaxa"
absworkdir = os.path.abspath(workdir)

configfi = "tests/data/test.config"

def test_0():
    if os.path.isfile("tests/data/precooked/otol_scraper.p"):
        # physcraper.debug(os.getcwd())
        conf = physcraper.ConfigObj(configfi, interactive=False)
        # physcraper.debug("conf")
        conf.unmapped = 'keep'
        # physcraper.debug("set unmapped")
        data_obj = pickle.load(open("tests/data/precooked/otol_tiny_dataobj.p", 'rb'))
        data_obj.workdir = absworkdir
        # physcraper.debug("dataobj loaded")
        ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
        ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/otol_tiny_gi_map.p", "rb"))
        # physcraper.debug("ids loaded")
        scraper = pickle.load(open("tests/data/precooked/otol_scraper.p", "rb"))
        # physcraper.debug("scraper loaded")
        # scraper2 = pickle.load(open("tests/data/precooked/otol_scraper.p", "rb"))
        num_keep = len(scraper.data.aln.taxon_namespace)
        # physcraper.debug('num_keep')

        # physcraper.debug(num_keep)
    # except:
    else:
        sys.stdout.write("\n\n No files present\n\n")
        conf = physcraper.ConfigObj(configfi)
        conf.unmapped = 'keep'
        aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
        data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                    config_obj=conf,
                                    study_id = study_id,
                                    tree_id = tree_id,
                                    phylesystem_loc = conf.phylesystem_loc)
        # physcraper.debug(len(data_obj.aln.taxon_namespace))
        pickle.dump(data_obj, open("tests/data/precooked/otol_tiny_dataobj.p", "wb" ))
        ids =  physcraper.IdDicts(conf, workdir=workdir)
        # physcraper.debug(os.getcwd())
        pickle.dump(ids.acc_ncbi_dict, open("tests/data/precooked/otol_tiny_gi_map.p", "wb"))
        data_obj.write_files()
        scraper = physcraper.PhyscraperScrape(data_obj, ids)
        # physcraper.debug(len(scraper.data.aln.taxon_namespace))
        # physcraper.debug("scraper obj made")
        pickle.dump(scraper.config, open("tests/data/precooked/otol_conf.p", "wb"))
        pickle.dump(scraper, open("tests/data/precooked/otol_scraper.p", "wb"))
        num_keep = len(scraper.data.aln.taxon_namespace)
        # physcraper.debug(num_keep)




    conf = physcraper.ConfigObj(configfi, interactive=False)
    conf.unmapped = 'remove'


    data_obj = pickle.load(open("tests/data/precooked/otol_tiny_dataobj.p", 'rb'))
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
