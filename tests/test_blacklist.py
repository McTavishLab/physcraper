import sys
import os
import json
import pickle
import shutil
from physcraper import wrappers, ConfigObj, IdDicts, FilterBlast, generate_ATT_from_files, AlignTreeTax
#


#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"



workdir="tests/output/test_blacklist"
configfi = "tests/data/test.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
# downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None

## make one run without blacklist
blacklist = None
noblack = os.path.join(workdir, "noblacklist")
absworkdir = os.path.abspath(noblack)

if not os.path.exists(os.path.join(absworkdir, "current_blast_run/")):
    os.makedirs(os.path.join(absworkdir, "current_blast_run/"))

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


noblackScrape =  FilterBlast(data_obj, ids)
noblackScrape._blasted = 1

src = "tests/data/precooked/fixed/tte_blast_files"

src_files = os.listdir(src)
for file_name in src_files:
    dest =  os.path.join(absworkdir, "current_blast_run/")
    print(dest)
    full_file_name = os.path.join(src, file_name)
    if (os.path.isfile(full_file_name)):
        shutil.copy(full_file_name, dest)


noblackScrape.read_blast()
noblackScrape.remove_identical_seqs()
        
noblackScrape.generate_streamed_alignment()


## one run with blacklist

blacklist = [429489230]


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
filteredScrape._blasted = 1
if not os.path.exists(os.path.join(absworkdir, "current_blast_run/")):
    os.makedirs(os.path.join(absworkdir, "current_blast_run/"))

src = "tests/data/precooked/fixed/tte_blast_files"

src_files = os.listdir(src)
for file_name in src_files:
    dest =  os.path.join(absworkdir, "current_blast_run/")
    # print(dest)
    full_file_name = os.path.join(src, file_name)
    if (os.path.isfile(full_file_name)):
        shutil.copy(full_file_name, dest)


filteredScrape.read_blast()
filteredScrape.remove_identical_seqs()
        
filteredScrape.generate_streamed_alignment()

for item in blacklist:
    for tax in filteredScrape.data.tre.taxon_namespace:
        # print(filteredScrape.data.otu_dict[tax.label])
        gi_id = filteredScrape.data.otu_dict[tax.label].get("^ncbi:gi")
        # print(tax, gi_id)

        try:
            assert item != gi_id
        except:
            print("test part1 failed")


    for tax in noblackScrape.data.tre.taxon_namespace:
        # print(filteredScrape.data.otu_dict[tax.label])
        gi_id = noblackScrape.data.otu_dict[tax.label].get("^ncbi:gi")
        # print(tax, gi_id)

        if item == gi_id:
            print("seq was not added in blacklist run")
            print("test works")

# test if it removes blacklist gi from already added aln:

for item in blacklist:
    for tax in noblackScrape.data.tre.taxon_namespace:
        # print(filteredScrape.data.otu_dict[tax.label])
        gi_id = noblackScrape.data.otu_dict[tax.label].get("^ncbi:gi")
        # print(tax, gi_id)

        if item == gi_id:
            print("blacklist gi was added in previous run")

print("now we want to remove it.")

blacklist = [429489230]

len_before = (len(noblackScrape.data.tre.taxon_namespace))
noblackScrape.blacklist = blacklist
noblackScrape.generate_streamed_alignment()
print(len_before)
print(len(noblackScrape.data.tre.taxon_namespace))

try:
    assert len_before-1 == len(noblackScrape.data.tre.taxon_namespace)
    print("pass test")
except:
    print("test failed")

