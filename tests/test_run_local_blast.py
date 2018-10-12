import os
import sys
import pickle
#from physcraper import FilterBlast, ConfigObj, IdDicts
import physcraper
import physcraper.local_blast as local_blast


sys.stdout.write("\ntests run_local_blast\n")


# tests if I can run a local blast query
workdir = "tests/output/test_run_local_blast"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)

try:
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape =  physcraper.FilterBlast(data_obj, ids)

blast_db = "otuSlagascanus"
blast_seq = "otuSlagascanus"

if not os.path.exists("{}/blast".format(filteredScrape.data.workdir)):
    os.makedirs("{}/blast/".format(filteredScrape.data.workdir))
path1 = '{}/tests/data/precooked/fixed/select-blast/*'.format(os.getcwd())

path2 = "{}/blast/".format(filteredScrape.data.workdir)
cmd = 'cp -r ' + path1 + ' ' + path2
os.system(cmd)

local_blast.run_local_blast(filteredScrape.data.workdir, blast_seq, blast_db)
blast_out = "{}/blast/output_otuSlagascanus_tobeblasted.xml".format(workdir)

if os.path.exists(blast_out):
    open(blast_out)
    sys.stdout.write("\ntest passed\n")
else:
    sys.stderr.write("\ntest failed\n")
