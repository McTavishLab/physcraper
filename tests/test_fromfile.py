
import sys
from physcraper import generate_ATT_from_files, generate_ATT_from_run


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
                                     ingroup_mrca=None)

    data_obj == True


def test_generate_ATT_from_run():
    workdir="tests/data/precooked/output"

    sys.stdout.write("\nTesting 'generate_ATT_from_run (fromfile.py)'\n")
    data_obj = generate_ATT_from_run(workdir=workdir)
      