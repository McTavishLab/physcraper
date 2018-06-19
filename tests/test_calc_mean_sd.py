import sys
import os
import json
import pickle
from math import sqrt
from Bio.Blast import NCBIXML
from physcraper import debug, wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast



print('test calculate_mean_sd')
# print("NOTE: the one which does reload the data, is not running, as it does not find a file....")

# ###########################
# print("this runs though:")
# print("It only runs, if the data are not reloaded, if reload it has the wrong workdir!")
# print('test calculate_mean_sd')

seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/mean_sd_test"
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
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)

    json.dump(otu_json, open(otu_jsonfi,"w"))

if os.path.isfile("{}/mean_sd_test.p".format(workdir)): 
    print("reload file setup")
    filteredScrape = pickle.load(open("{}/mean_sd_test.p".format(workdir),'rb'))
 
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
    filteredScrape.dump('{}/mean_sd_test.p'.format(workdir))




# test begins
print("test begins")
fn = 'Senecio_scopolii_subsp._scopolii'

# partly copy of read_local_blast
general_wd = os.getcwd()

os.chdir(os.path.join(filteredScrape.workdir, "blast"))
debug(os.getcwd())
output_blast = "output_{}_tobeblasted.xml".format(fn)
xml_file = open(output_blast)
os.chdir(general_wd)
blast_out = NCBIXML.parse(xml_file)
hsp_scores = {}
                
add_hsp = 0
for record in blast_out:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            # filter by e-value
            ## !!! maybe don"t do that....
            if hsp.expect < 0.003:
                gi = int(alignment.title.split(" ")[1])
                hsp_scores[gi] = {"hsp.bits" : hsp.bits, "hsp.score" : hsp.score, "alignment.length" : alignment.length, "hsp.expect" : hsp.expect}
                add_hsp = add_hsp + float(hsp.bits)
                #debug(add_hsp)
# make values to select for blast search, calculate standard deviation, mean
mean_sed = filteredScrape.calculate_mean_sd(hsp_scores)


sum_hsp = len(hsp_scores)
mean = (add_hsp/sum_hsp)

sd_all = 0
for item in hsp_scores:
    #debug(item)
    val = hsp_scores[item]["hsp.bits"]
    sd = (val-mean)*(val-mean)
    sd_all += sd
sd_val = sqrt(sd_all/sum_hsp)


try: 
    assert round(sd_val, 4) == round(mean_sed['sd'], 4)
    assert round(mean,4) == round(mean_sed['mean'], 4)
    print('mean and sd are the same')
except:
    print('mean or sd are different')


# this does not run:
# ###########################

# seqaln = "tests/data/tiny_test_example/test.fas"
# mattype = "fasta"
# trfn = "tests/data/tiny_test_example/test.tre"
# schema_trf = "newick"
# workdir = "tests/output/mean_test"
# configfi = "tests/data/blubb_localblast.config"
# id_to_spn = r"/home/blubb/Documents/gitdata/physcraper/tests/data/tiny_test_example/test_nicespl.csv"
# otu_jsonfi = "{}/otu_dict.json".format(workdir)
# treshold = 2
# selectby = 'blast'
# downtorank = None
# add_local_seq = None
# id_to_spn_addseq_json = None


# try:
#     data_obj = pickle.load(open("tests/data/tiny_dataobj.p", 'rb'))
#     conf = ConfigObj(configfi)
#     ids = IdDicts(conf, workdir=data_obj.workdir)
#     ids.gi_ncbi_dict = pickle.load(open("tests/data/tiny_gi_map.p", "rb" ))
#     filteredScrape =  FilterBlast(data_obj, ids)
#     filteredScrape._blasted = 1
#     filteredScrape.dump('{}/mean_sd_test.p'.format(workdir))




# # test begins

#     fn = 'Senecio_scopolii_subsp._scopolii'

#     # partly copy of read_local_blast
#     general_wd = os.getcwd()
#     os.chdir(os.path.join(filteredScrape.workdir, "blast"))
#     output_blast = "output_{}_tobeblasted.xml".format(fn)
#     xml_file = open(output_blast)
#     os.chdir(general_wd)
#     blast_out = NCBIXML.parse(xml_file)
#     hsp_scores = {}
                    
#     add_hsp = 0
#     for record in blast_out:
#         for alignment in record.alignments:
#             for hsp in alignment.hsps:
#                 # filter by e-value
#                 ## !!! maybe don"t do that....
#                 if hsp.expect < 0.003:
#                     gi = int(alignment.title.split(" ")[1])
#                     hsp_scores[gi] = {"hsp.bits" : hsp.bits, "hsp.score" : hsp.score, "alignment.length" : alignment.length, "hsp.expect" : hsp.expect}
#                     add_hsp = add_hsp + float(hsp.bits)
#                     #debug(add_hsp)
#     # make values to select for blast search, calculate standard deviation, mean
#     mean_sed = filteredScrape.calculate_mean_sd(hsp_scores)


#     sum_hsp = len(hsp_scores)
#     mean = (add_hsp/sum_hsp)

#     sd_all = 0
#     for item in hsp_scores:
#         #debug(item)
#         val = hsp_scores[item]["hsp.bits"]
#         sd = (val-mean)*(val-mean)
#         sd_all += sd
#     sd_val = sqrt(sd_all/sum_hsp)


#     try: 
#         assert round(sd_val, 4) == round(mean_sed['sd'], 4)
#         assert round(mean,4) == round(mean_sed['mean'], 4)
#         print('mean and sd are the same')
#     except:
#         print('mean or sd are different')

# except:
#     if not os.path.exists("{}".format(workdir)):
#         os.makedirs("{}".format(workdir))


#     otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)

#     with open(otu_jsonfi,"w") as outfile:
#         json.dump(otu_json, outfile)


#     conf = ConfigObj(configfi)
#     data_obj = generate_ATT_from_files(seqaln=seqaln, 
#                                  mattype=mattype, 
#                                  workdir=workdir,
#                                  treefile=trfn,
#                                  schema_trf = schema_trf,
#                                  otu_json=otu_jsonfi,
# #                                 email = conf.email,
#                                  ingroup_mrca=None)

#     filteredScrape =  FilterBlast(data_obj, ids)
#     filteredScrape._blasted = 1
    
#     filteredScrape.dump('{}/mean_sd_test.p'.format(workdir))
#     print('rerun the test, the files were not present')


