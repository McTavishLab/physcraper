import sys
import os
import json
from physcraper import generate_ATT_from_files, AlignTreeTax, OtuJsonDict, ConfigObj, IdDicts, PhyscraperScrape
from pytest import mark
from physcraper.opentree_helpers import bulk_tnrs_load
from physcraper.helpers import cd
from dendropy import DnaCharacterMatrix
from physcraper import ncbi_data_parser



web = mark.web

def test_cd():
    workdir="tests/output/"
    with cd(workdir):
        assert(os.getcwd() == os.path.abspath(workdir))





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
    alnfile =  "tests/data/tiny_test_example/test.fas"
    aln_schema = "fasta"
    otu_json = "tests/data/tiny_test_example/main.json"
    workdir = "tests/data/precooked/tiny_local"
    data_obj = generate_ATT_from_files(workdir= workdir,
                                        configfile=None,
                                        alnfile = alnfile,
                                        aln_schema = aln_schema,
                                        treefile = treefile,
                                        otu_json = otu_json,
                                        tree_schema = tree_schema)
    ids = IdDicts()
    ids.ncbi_parser = ncbi_data_parser.Parser(names_file="taxonomy/names.dmp",
                                               nodes_file="taxonomy/nodes.dmp")
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.config.blast_loc = "local"
    scraper._blasted = 1
    scraper.read_blast_wrapper(blast_dir="tests/data/precooked/tiny_local/blast_run_test")
    assert(len(scraper.new_seqs_otu_id)==17)
