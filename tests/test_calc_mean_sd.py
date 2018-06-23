import sys
import os
import json
import pickle
from math import sqrt
from Bio.Blast import NCBIXML
from physcraper import debug, wrappers, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast



print('test calculate_mean_sd')


seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/mean_sd_test"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = 'blast'
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None




absworkdir = os.path.abspath(workdir)


try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)




# test begins
print("test begins")
fn = 'Senecio_scopolii_subsp._scopolii'

# partly copy of read_local_blast
general_wd = os.getcwd()
if not os.path.exists(os.path.join(filteredScrape.workdir, "blast")):
    os.makedirs(os.path.join(filteredScrape.workdir, "blast"))
os.chdir(os.path.join(filteredScrape.workdir, "blast"))
debug(os.getcwd())

fn_path = '/home/blubb/Documents/gitdata/physcraper/tests/data/precooked/fixed/local-blast/{}'.format(fn)
filteredScrape.run_local_blast(fn_path, fn_path, output=os.path.join(filteredScrape.workdir, "blast/output_{}.xml".format(fn)) )

output_blast = os.path.join(filteredScrape.workdir, "blast/output_{}.xml".format(fn))
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
            # if hsp.expect < 0.003:
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

print((sd_val, 4), round(mean_sed['sd'], 4))
print(mean,4), round(mean_sed['mean'], 4)
try: 
    assert round(sd_val, 4) == round(mean_sed['sd'], 4)
    assert round(mean,4) == round(mean_sed['mean'], 4)
    print('mean and sd are the same')
except:
    print('mean or sd are different')

