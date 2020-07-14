import sys
import os
import json
from physcraper import generate_ATT_from_files, AlignTreeTax, OtuJsonDict, ConfigObj, IdDicts
from pytest import mark
from physcraper.opentree_helpers import bulk_tnrs_load
from physcraper.helpers import cd


web = mark.web

def test_cd():
    workdir="tests/output/"
    with cd(workdir):
        assert(os.getcwd() == os.path.abspath(workdir))

def test_owndata():
    seqaln= "tests/data/tiny_test_example/test.fas"
    mattype="fasta"
    trfn= "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    workdir="tests/output/owndata"
    configfi = "tests/data/test.config"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    """Tests if your own input files will generate a data object of class AlignTreeTax
    """

    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    ids = IdDicts()
    if os.path.exists(otu_jsonfi):
        print("load json")
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = OtuJsonDict(id_to_spn, ids)
        json.dump(otu_json, open(otu_jsonfi,"w"))

    data_obj = generate_ATT_from_files(alnfile=seqaln,
                                 aln_schema=mattype,
                                 workdir=workdir,
                                 configfile=configfi,
                                 treefile=trfn,
                                 tree_schema = schema_trf,
                                 otu_json=otu_jsonfi,
                                 search_taxon=None)


    assert isinstance(data_obj, AlignTreeTax)


import physcraper
from dendropy import DnaCharacterMatrix


def test_load_json():
    inputfi = "tests/data/bulk_tnrs.json"
    otu_dict = bulk_tnrs_load(inputfi)
    assert len(otu_dict) == 16
    print(otu_dict)



def test_owndata_bulktnrs():
    seqaln= "tests/data/tiny_test_example/test.fas"
    mattype="fasta"
    trfn= "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    workdir="tests/output/owndata"
    configfi = "tests/data/test.config"
    otu_jsonfi = "tests/data/tiny_test_example/main.json"

    """Tests if your own input files will generate a data object of class AlignTreeTax
    """

    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)

    otu_dict = bulk_tnrs_load(otu_jsonfi)  
    print(otu_dict)
    data_obj = generate_ATT_from_files(alnfile=seqaln,
                                 aln_schema=mattype,
                                 workdir=workdir,
                                 configfile=configfi,
                                 treefile=trfn,
                                 tree_schema = schema_trf,
                                 otu_json=otu_dict,
                                 search_taxon=None)


    assert isinstance(data_obj, AlignTreeTax)


test_owndata_bulktnrs()


def test_add_all_local():
    treefile = "tests/data/tiny_test_example/test.tre"
    tree_schema = "newick"
    alnfile =  "tests/data/tiny_test_example/oneseq.fas"
    aln_schema = "fasta"
     otu_json = "tests/data/tiny_test_example/main.json"
    workdir = "tests/data/precooked/tiny_local"
    data_obj = generate_ATT_from_files(workdir= workdir,
                                        configfile=conf,
                                        alnfile = alnfile,
                                        aln_schema = aln_schema,
                                        treefile = treefile,
                                        otu_json = otu_dict,
                                        tree_schema = args.tree_schema,
                                        search_taxon=search_ott_id)
    ids = physcraper.IdDicts()
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    scraper.blast_loc = "local"
    scraper.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tiny_local/blast_run_test")
    scraper.remove_identical_seqs()
    assert(len(scraper.new_seqs)==17)
