from physcraper import AlignTreeTax, generate_ATT_from_files, ConfigObj, IdDicts,  PhyscraperScrape
import sys
import os
from physcraper import wrappers_nonstandard
#
seqaln= "docs/owndata/senecio_its.fasta"
mattype="fasta"
trfn= "docs/owndata/its_new.tre"
schema_trf = "newick"
workdir="owndata_testrun"
otujson = "docs/owndata/ott_info_owntree_its.txt"
configfi = "example.config"
id_to_spn = "docs/owndata/uniquetip_to_name_its.csv"
cwd = os.getcwd()  



"""Tests if your own input files will generate a data object of class AlignTreeTax
"""

sys.stdout.write("setting up Data Object\n")
sys.stdout.flush()
#read the config file into a configuration object
conf = ConfigObj(configfi)




otu_json = wrappers_nonstandard.OtuJsonDict(id_to_spn, configfi)




#Generate an linked Alignment-Tree-Taxa object
data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                 mattype=mattype, 
                                 workdir=workdir,
                                 treefile=trfn,
                                 schema_trf = schema_trf,
                                 otu_json=otu_json,
                                 ingroup_mrca=None)

print(type(data_obj))
if isinstance(data_obj, AlignTreeTax):
    print data_obj, "is of type AlignTreeTax. Success, your input should be working."