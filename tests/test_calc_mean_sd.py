import sys
import os
import pickle
from math import sqrt
from Bio.Blast import NCBIXML
from physcraper import ConfigObj, IdDicts, FilterBlast

sys.stdout.write("\ntests calculate_mean_sd\n")

workdir = "tests/output/mean_sd_test"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)

try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape = FilterBlast(data_obj, ids)

# test begins
fn = 'Senecio_scopolii_subsp._scopolii'
# partly copy of read_local_blast
general_wd = os.getcwd()
if not os.path.exists(os.path.join(filteredScrape.workdir, "blast")):
    os.makedirs(os.path.join(filteredScrape.workdir, "blast"))
os.chdir(os.path.join(filteredScrape.workdir, "blast"))
fn_path = '/home/blubb/Documents/gitdata/physcraper/tests/data/precooked/fixed/local-blast/{}'.format(fn)
filteredScrape.run_local_blast(fn_path, fn_path,
                               output=os.path.join(filteredScrape.workdir, "blast/output_{}.xml".format(fn)))
output_blast = os.path.join(filteredScrape.workdir, "blast/output_{}.xml".format(fn))
xml_file = open(output_blast)
os.chdir(general_wd)
blast_out = NCBIXML.parse(xml_file)
hsp_scores = {}            
add_hsp = 0
for record in blast_out:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            gi = int(alignment.title.split(" ")[1])
            hsp_scores[gi] = {"hsp.bits": hsp.bits, "hsp.score": hsp.score, "alignment.length": alignment.length, "hsp.expect": hsp.expect}
            add_hsp = add_hsp + float(hsp.bits)
# make values to select for blast search, calculate standard deviation, mean
mean_sed = filteredScrape.calculate_mean_sd(hsp_scores)
sum_hsp = len(hsp_scores)
mean = (add_hsp / sum_hsp)
sd_all = 0
for item in hsp_scores:
    val = hsp_scores[item]["hsp.bits"]
    sd = (val - mean) * (val - mean)
    sd_all += sd
sd_val = sqrt(sd_all / sum_hsp)
# print((sd_val, 4), round(mean_sed['sd'], 4))
# print(mean,4), round(mean_sed['mean'], 4)
try:
    assert round(sd_val, 4) == round(mean_sed['sd'], 4)
    assert round(mean, 4) == round(mean_sed['mean'], 4)
    sys.stdout.write('\nTEST PASSED: mean and sd are the same\n')
except:
    # print('TEST FAILED: mean or sd are different')
    sys.stderr.write("\ntest failed\n")

