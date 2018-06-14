import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast
import pickle#

# tests if the building of sp_d and sp_seq_dict is working correclty
#
seqaln = "tests/data/non_unicode_names/test.fas"
mattype = "fasta"
trfn = "tests/data/non_unicode_names/test.tre"
schema_trf = "newick"
workdir = "tests/output/edit_dict_key"
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/non_unicode_names/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None

print("testing: edit_dict_key")

try:
    if os.path.exists(otu_jsonfi):
        print("reload from file")
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        json.dump(otu_json, open(otu_jsonfi,"w"))

    conf = ConfigObj(configfi)
           
    # #aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

    #Generate an linked Alignment-Tree-Taxa object
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                     mattype=mattype, 
                                     workdir=workdir,
                                     treefile=trfn,
                                     schema_trf = schema_trf,
                                     otu_json=otu_jsonfi,
    #                                 email = conf.email,
                                     ingroup_mrca=None)


    ids = IdDicts(conf, workdir=workdir)
    #Now combine the data, the ids, and the configuration into a single physcraper scrape object
    filteredScrape =  FilterBlast(data_obj, ids)


    print("run test")
  
    filteredScrape.data.otu_dict = {str(k): str(v) for k, v in filteredScrape.data.otu_dict.items()}

    tax_namespace1 = filteredScrape.data.otu_dict.keys()
    filteredScrape.data.edit_dict_key()
    tax_namespace2 = filteredScrape.data.otu_dict.keys()
    assert sorted(tax_namespace1) != sorted(tax_namespace2)
    print("edit_dict_key changes the names")

except:
    print("Test edit_dict_key is not running as expected")

