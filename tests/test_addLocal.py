import pickle  #
from physcraper.configobj import ConfigObj
from physcraper.ids import IdDicts
from physcraper.opentree_helpers import OtuJsonDict
from physcraper.filterblast import FilterBlast
import os
import json
import sys
import pickle
from pytest import mark

localblast = mark.localblast



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

threshold = 10
selectby = "blast"
downto = None
add_local_seq = "tests/data/local_seqs"
id_to_spn_addseq = "tests/data/tipnTOspn_localAdd.csv"

###################
absworkdir = os.path.abspath(workdir)
conf = ConfigObj("tests/data/test.config", interactive=False)


@localblast
def test_add_local():
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.blacklist = blacklist

    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    if os.path.exists(otu_jsonfi_local):
        otu_json_local = json.load(open(otu_jsonfi_local))
    else:
        otu_json_local = OtuJsonDict(id_to_spn_addseq, ids)
        json.dump(otu_json_local, open(otu_jsonfi_local, "w"))

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
        filteredScrape.run_blast_wrapper()
        filteredScrape.read_blast_wrapper()
        filteredScrape.remove_identical_seqs()

    test = False
    for key in filteredScrape.data.otu_dict.keys():
        if '^ncbi:title' in filteredScrape.data.otu_dict[key].keys():
            if filteredScrape.data.otu_dict[key]['^ncbi:title'] == "unpublished":
                test = True
                break
    assert test == True

