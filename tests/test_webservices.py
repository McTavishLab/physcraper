import subprocess
import os
import json

from physcraper import OtuJsonDict, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from physcraper import treetaxon, opentree_helpers
#


seqaln= "tests/data/tiny_test_example/oneseq.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"

workdir="tests/data/tmp/webservices"

configfi = "tests/data/test.config"
id_to_spn = "tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

"""Generates the files needed for the tests.
"""

if not os.path.exists("{}".format(workdir)):
    os.makedirs("{}".format(workdir))

conf_base = ConfigObj(configfi)

ids_base = IdDicts(configfi)

otu_json = OtuJsonDict(id_to_spn, ids_base)
with open(otu_jsonfi,"w") as outfile:
    json.dump(otu_json, outfile)

"""Extract ott ids from `otu_json` dictionary (???):
"""

ottids = [otu_json[ite]['^ot:ottId'] for ite in otu_json]


"""To get the MRCA of the matched taxa, use function `opentree_helpers.get_mrca_ott()`:
"""

mrca = opentree_helpers.get_mrca_ott(ottids)


"""Generate the alignment-taxon object:
"""

data_obj_base = generate_ATT_from_files(alnfile=seqaln,
                             aln_schema=mattype,
                             workdir=workdir,
                             configfile=configfi,
                             treefile=trfn,
                             tree_schema = schema_trf,
                             otu_json=otu_jsonfi,
                             search_taxon=mrca)


scraper = PhyscraperScrape(data_obj_base, ids_base)

gb_id = 'JX895264.1'
def test_get_full_seq():
    taxid,taxname, seq = ids_base.get_tax_seq_acc(gb_id)
    assert(taxid == '1268581')
    assert(len(seq) == 752)
    match = 'cagcattgttccaaa'
    revseq = scraper.check_complement(match, seq, gb_id)
    assert(revseq.startswith('ccttcattttcagcattgttccaa'))


def test_get_from_treebase():
    subprocess.check_call(["python", "bin/physcraper_run.py", 
                            "-s" "pg_55",
                            "-t", "tree5864",
                            "-tb",
                            "-no_est",
                            "-o", "tests/data/tmp/pg_55_treebase"])

def test_web_blast():
    scraper.run_blast_wrapper()
    scraper.read_blast_wrapper()
    assert(scraper)


def test_reroot():
    tr = treetaxon.generate_TreeTax_from_run(workdir="tests/data/precooked/ot_350")
    opentree_helpers.root_tree_from_synth(tree=tr.tre, otu_dict=tr.otu_dict)


def test_find_trees():
    subprocess.check_call(["python", "bin/find_trees.py",
                            "--taxon_name", "Malvaceae",
                            "--treebase",
                            "-o", "tests/tmp/malvacea.txt"])
    assert(os.path.exists("tests/tmp/malvacea.txt"))
    subprocess.check_call(["python", "bin/find_trees.py",
                            "--ott_id", "124219",
                            "--treebase",
                            "-o", "tests/tmp/orcinus.txt"])
    assert(os.path.exists("tests/tmp/orcinus.txt"))



