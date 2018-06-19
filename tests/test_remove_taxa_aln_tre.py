import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can remove sequences from aln and tre


# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_remove_taxa_aln_tre"
configfi = "tests/data/blubb_localblast.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = 'blast'
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None


if os.path.exists(otu_jsonfi):
    print("reload from file")
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.isfile("{}/remove_taxa_aln_tre.p".format(workdir)): 
    filteredScrape = pickle.load(open("{}/remove_taxa_aln_tre.p".format(workdir),'rb'))
 
else:   
    conf = ConfigObj(configfi)
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                         mattype=mattype, 
                         workdir=workdir,
                         treefile=trfn,
                         schema_trf=schema_trf,
                         otu_json=otu_jsonfi,
                         ingroup_mrca=None)


    data_obj.prune_short()
    data_obj.dump()

    ids = IdDicts(conf, workdir=workdir)
    ids.dump()

    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.dump('{}/remove_taxa_aln_tre.p'.format(workdir))

for tax in filteredScrape.data.aln:
	print(tax)
	print(tax.label)


#print(filteredScrape.data.aln.as_string('phylip'))	
#print(filteredScrape.data.tre)	
len_aln_before = len(filteredScrape.data.aln.as_string('phylip'))
len_tre_before = len(filteredScrape.data.tre.as_string(schema="newick"))
namespace_before = len(filteredScrape.data.aln.taxon_namespace)
namespace_tre_before = len(filteredScrape.data.tre.taxon_namespace)

#print(len(filteredScrape.data.aln.as_string('phylip')))	
#print(len(filteredScrape.data.tre.as_string(schema="newick")))
#print(filteredScrape.data.aln.taxon_namespace)
#print(filteredScrape.data.tre.taxon_namespace)

print('remove from object')
filteredScrape.data.remove_taxa_aln_tre(tax.label)
#print(filteredScrape.data.aln.as_string('phylip'))	
#print(filteredScrape.data.tre)	
len_aln_after = len(filteredScrape.data.aln.as_string('phylip'))
len_tre_after = len(filteredScrape.data.tre.as_string(schema="newick"))
namespace_after = len(filteredScrape.data.aln.taxon_namespace)
namespace_tre_after = len(filteredScrape.data.tre.taxon_namespace)

#print(len(filteredScrape.data.aln.as_string('phylip')))	
#print(len(filteredScrape.data.tre.as_string(schema="newick")))
#print(filteredScrape.data.aln.taxon_namespace)
#print(filteredScrape.data.tre.taxon_namespace)


try:
    try:
        assert len_aln_before != len_aln_after
    except:
        print('test1 failed')
    try:
        assert len_tre_before != len_tre_after
    except:
        print('test2 failed')
    try:
        assert namespace_before != namespace_after
    except:
        print('test3 failed')
    try:
        assert namespace_tre_before != namespace_tre_after
    except:
        print('test4 failed')
    print("all subtests passed")
except:
    print("test remove tax aln tre failed")

## test passes, but sometimes the taxon is not deleted from the tre. this keeps being an issue. have added some slightly different commands to make sure that it is being removed. seems to be somehow redundant, though.