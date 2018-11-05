from physcraper import OtuJsonDict, ConfigObj, IdDicts, FilterBlast
import os
import json
import sys
import pickle

#################################

workdir = "tests/output/test_addLocal"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)

seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
blacklist = None

id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
otu_jsonfi_local = "{}/otu_dict_local.json".format(workdir)

cwd = os.getcwd()
threshold = 10
selectby = "blast"
downto = None
add_local_seq = "tests/data/local_seqs"
id_to_spn_addseq = "tests/data/tipnTOspn_localAdd.csv"

###################

try:
    conf = ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest setup failed\n\n")
    sys.exit(-1)

if not os.path.exists("{}".format(workdir)):
    os.makedirs("{}".format(workdir))

if os.path.exists(otu_jsonfi_local):
    otu_json_local = json.load(open(otu_jsonfi_local))
else:
    otu_json_local = OtuJsonDict(id_to_spn_addseq, ids)
    json.dump(otu_json_local, open(otu_jsonfi_local, "w"))


sys.stdout.write("\ntest addLocal\n")

# Prune sequences below a certain length threshold
data_obj.prune_short()
data_obj.write_files()
data_obj.write_labelled(label='^ot:ottTaxonName', add_gb_id=True)
data_obj.write_otus("otu_info", schema='table')
data_obj.dump()

sys.stdout.write("setting up id dictionaries\n")
sys.stdout.flush()

ids = IdDicts(conf, workdir=workdir)

# Now combine the data, the ids, and the configuration into a single physcraper scrape object
filteredScrape = FilterBlast(data_obj, ids)
filteredScrape.blacklist = blacklist


if add_local_seq is not None:
    filteredScrape.unpublished = True
if filteredScrape.unpublished is True:  # use unpublished data
    # filteredScrape.unpublished = True
    filteredScrape.data.unpubl_otu_json = otu_json_local
    filteredScrape.write_unpubl_blastdb(add_local_seq)

    # filteredScrape.make_otu_dict_entry_unpubl()
    filteredScrape.run_blast()
    filteredScrape.read_blast_wrapper()
    filteredScrape.remove_identical_seqs()

test = False
for key in filteredScrape.data.otu_dict.keys():
    if '^ncbi:title' in filteredScrape.data.otu_dict[key].keys():
        if filteredScrape.data.otu_dict[key]['^ncbi:title'] == "unpublished":
            test = True
            break
if test:
    sys.stdout.write("\ntest passed\n")
else:
    # print("test remove tax aln tre failed")
    sys.stderr.write("\ntest failed\n")
