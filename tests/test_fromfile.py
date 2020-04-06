
import json
import sys
from physcraper import generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from dendropy import Tree,\
                     DnaCharacterMatrix,\
                     DataSet,\
                     datamodel

#Use OpenTree phylesystem identifiers to get study and tree

def test_generate_ATT_from_file():

    seqaln = "tests/data/input.fas"
    mattype="fasta"
    workdir="tests/fromfile"
    treefile = "tests/data/input.tre"
    otu_jsonfi = "tests/data/otu_dict.json"
    schema_trf = "newick"
    configfi = "tests/data/test.config"

    sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")

    conf = ConfigObj(configfi, interactive=False)

    data_obj = generate_ATT_from_files(alnfile=seqaln, 
                                     aln_schema=mattype, 
                                     workdir=workdir,
                                     configfile=configfi,
                                     treefile=treefile,
                                     tree_schema=schema_trf,
                                     otu_json=otu_jsonfi,
                                     ingroup_mrca=None)



    data_obj == True
      