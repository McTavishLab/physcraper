
import sys
from physcraper import IdDicts, PhyscraperScrape, generate_ATT_from_files, generate_ATT_from_run


def test_generate_ATT_from_files():

    seqaln = "tests/data/input.fas"
    mattype="fasta"
    workdir="tests/fromfile"
    treefile = "tests/data/input.tre"
    otu_jsonfi = "tests/data/otu_dict.json"
    schema_trf = "newick"
    configfi = "tests/data/test.config"

    sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")
    data_obj = generate_ATT_from_files(alnfile=seqaln, 
                                     aln_schema=mattype, 
                                     workdir=workdir,
                                     configfile=configfi,
                                     treefile=treefile,
                                     tree_schema=schema_trf,
                                     otu_json=otu_jsonfi,
                                     search_taxon=None)

    data_obj == True


def test_generate_ATT_from_run():
    workdir="tests/data/precooked/output"

    sys.stdout.write("\nTesting 'generate_ATT_from_run '\n")
    data_obj = generate_ATT_from_run(workdir=workdir)
      

def test_example():
    indir="docs/examples/pg_55_web"
    workdir = "tests/tmp/example_test"
    data_obj = generate_ATT_from_run(workdir=indir, start_files = "input")
    data_obj.workdir = workdir
    ids = IdDicts(data_obj.config)
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.read_blast_wrapper("docs/examples/pg_55_web/blast_run_pg_55tree5864/")
    assert(len(scraper.new_seqs_otu_id) == 25)