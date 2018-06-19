
import json
import sys
from physcraper import generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from dendropy import Tree,\
                     DnaCharacterMatrix,\
                     DataSet,\
                     datamodel

#Use OpenTree phylesystem identifiers to get study and tree
seqaln = "tests/data/input.fas"
mattype="fasta"
workdir="tests/fromfile"
treefile = "tests/data/input.tre"
otu_jsonfi = "tests/data/otu_dict.json"
schema_trf = "newick"

sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")
data_obj = generate_ATT_from_files(seqaln = seqaln,
                                   mattype = mattype,
                                       workdir =workdir,
                                       treefile = treefile,
                                       schema_trf = schema_trf,
                                       otu_json = otu_jsonfi)
if data_obj:
    sys.stdout.write("\nTest fromfile.py passed\n")
else:
    sys.stdout.write("\nTest fromfile.py FAILED'\n")
