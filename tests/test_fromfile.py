
import sys
from Bio.Blast import NCBIXML
from physcraper import IdDicts, PhyscraperScrape, ConfigObj, generate_ATT_from_files, generate_ATT_from_run, ncbi_data_parser


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
    assert(len(scraper.new_seqs_otu_id) == 3)


"""
def test_webbased_blast():
    indir="docs/examples/pg_55_web"
    workdir = "tests/tmp/example_test"
    data_obj = generate_ATT_from_run(workdir=indir, start_files = "input")
    data_obj.workdir = workdir
    ids = IdDicts(data_obj.config)
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.read_webbased_blast_query('docs/examples/pg_55_web/blast_run_pg_55tree5864/otu376425.xml')
    len(scraper.new_seqs)
    scraper.read_blast_wrapper('docs/examples/pg_55_web/blast_run_pg_55tree5864')

def test_local_blast():
    indir="docs/examples/pg_55_web"
    workdir = "tests/tmp/example_loc_test"
    data_objloc = generate_ATT_from_run(workdir=indir, start_files = "input")
    data_objloc.workdir = workdir
    data_objloc.config.blastdb = "/home/ejmctavish/ncbi/localblastdb"
    data_objloc.config.set_local()
    idsloc = IdDicts(data_obj.config)
    idsloc.ncbi_parser = ncbi_data_parser.Parser(names_file=data_objloc.config.ncbi_names,
                            nodes_file=data_objloc.config.ncbi_nodes)
    scraperloc = PhyscraperScrape(data_objloc, idsloc)
    scraperloc.read_local_blast_query('docs/examples/pg_55_local/blast_run_pg_55tree5864_ndhf/otu376425.txt')
    scraperloc.read_blast_wrapper('docs/examples/pg_55_local/blast_run_pg_55tree5864_ndhf')

len(scraperloc.new_seqs)

"""



"""
def test_check_complement():
indir="docs/examples/pg_55_web"
workdir = "tests/tmp/example_test"
data_obj = generate_ATT_from_run(workdir=indir, start_files = "input")
data_obj.workdir = workdir
ids = IdDicts(data_obj.config)
scraper = PhyscraperScrape(data_obj, ids)
scraper.data.gb_dict = {}
result_handle = open('docs/examples/pg_55_web/blast_run_pg_55tree5864/otu376425.xml')
blast_records = NCBIXML.parse(result_handle)
i = 0
ok = 0
for blast_record in blast_records:
    i += 1
    print(i)
    for alignment in blast_record.alignments:
            ok += 1
            if ok > 1:
                break
            for hsp in alignment.hsps:
                    print(hsp)
                    if float(hsp.expect) < float(scraper.config.e_value_thresh):
                        gb_id = alignment.title.split("|")[3]  # 1 is for gi
                        print(gb_id)
                        if len(gb_id.split(".")) == 1:
                             sys.stdout.write("skipping acc {}, incorrect format\n".format(gb_id))
                        elif gb_id not in scraper.data.gb_dict:  # skip ones we already have
                            taxid,taxname, seq = scraper.ids.get_tax_seq_acc(gb_id)
                            print(gb_id)
                            print(taxid)
                            print(seq)
                            gi_id = alignment.title.split('|')[1]
                            gb_acc = alignment.accession
                            stitle = alignment.title
                            print(match)
                            match = hsp.sbjct.lower().replace('-','')
                            if match not in seq.lower():
                                print("was mismatch")
                            else: 
                                    print("was match")"""

