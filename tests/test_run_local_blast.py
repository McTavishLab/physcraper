import os
import json
import pickle
from physcraper import FilterBlast, wrappers, ConfigObj, generate_ATT_from_files, IdDicts

# tests if I can run a local blast query


# I need to generate a FilterBlast object first
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_run_local_blast"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = None
selectby = None
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

# add var which you need to use for the test:
blast_db = "S_lagascanus"
blast_seq = "S_lagascanus"


if not os.path.exists("{}/blast".format(filteredScrape.data.workdir)):
    os.makedirs("{}/blast/".format(filteredScrape.data.workdir))
path1 = '/home/blubb/Documents/gitdata/physcraper/tests/data/precooked/fixed/'
path2 = "{}/blast/".format(filteredScrape.data.workdir)
os.system('cp -r' + path1 + ' ' + path2)

filteredScrape.run_local_blast(blast_seq, blast_db)

blast_out = "{}/blast/output_S_lagascanus_tobeblasted.xml".format(workdir)
# print(blast_out)
# xml_file = open(blast_out)
# print(xml_file)
if os.path.exists(blast_out):
    xml_file = open(blast_out)
    # print(xml_file)
    print("succesfully read local blast output file")
else:
    print("error")