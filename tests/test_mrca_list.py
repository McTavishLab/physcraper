# package import
import os
import json
import pickle
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts, FilterBlast

# define here your files
def test_mrca_list():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = [723076, 710505, 187044, 4727685, 4728090, 711399 ]

    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    conf.blast_loc='remote' #saves time over loading names and nodes, and they aren't used here
    ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)

    # print(ids.mrca_ott, ids.mrca_ncbi)

    assert len(ids.mrca_ncbi) >= 2
    assert ids.mrca_ott == ingroup_mrca
    assert ids.mrca_ott != ids.mrca_ncbi

    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids)
    assert len(filteredScrape.ids.mrca_ncbi) >= 2
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) == 38


def test_no_mrca():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = None
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    conf.blast_loc='remote' #saves time over loading names and nodes, and they aren't used here
    ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)

    # print(ids.mrca_ott, ids.mrca_ncbi)

    assert len(ids.mrca_ncbi) == 1
    assert ids.mrca_ott == ingroup_mrca
    assert ids.mrca_ott != ids.mrca_ncbi

    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids)
    assert len(filteredScrape.ids.mrca_ncbi) >= 2
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) <= 38

def test_higher_mrca():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/test_mrcalist_local"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)

    ingroup_mrca = 557768
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    conf.blast_loc='remote' #saves time over loading names and nodes, and they aren't used here
    ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)

    # print(ids.mrca_ott, ids.mrca_ncbi)

    assert len(ids.mrca_ncbi) >= 2
    assert ids.mrca_ott == ingroup_mrca
    assert ids.mrca_ott != ids.mrca_ncbi

    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    filteredScrape = FilterBlast(data_obj, ids)
    assert len(filteredScrape.ids.mrca_ncbi) >= 2
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) > 38


# #EJM version
# def test_mrca_list():
#     seqaln = "tests/data/tiny_test_example/test.fas"
#     mattype = "fasta"
#     trfn = "tests/data/tiny_test_example/test.tre"
#     schema_trf = "newick"
#     id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
#     workdir = "tests/output/test_mrcalist_local"
#     configfi = "tests/data/test.config"
#     otu_jsonfi = "{}/otu_dict.json".format(workdir)

#     ingroup_mrca = [723076, 710505, 187044, 4727685, 4728090, 711399]

#     # setup the run
#     if not os.path.exists("{}".format(workdir)):
#         os.makedirs("{}".format(workdir))

#     conf = ConfigObj(configfi)
#     conf.blast_loc='remote' #saves time over loading names and nodes, and they aren't used here
#     ids = IdDicts(conf, workdir=workdir)

#     # print(ids.mrca_ott, ids.mrca_ncbi)

#     # assert len(ids.mrca_ncbi) >= 2
#     # assert ids.mrca_ott == ingroup_mrca
#     # assert ids.mrca_ott != ids.mrca_ncbi
#     wrappers.make_otujsondict(id_to_spn, workdir, ids, local=False)
#     data_obj = wrappers.generate_ATT_from_files(seqaln=seqaln,
#                                            mattype=mattype,
#                                            workdir=workdir,
#                                            config_obj=conf,
#                                            treefile=trfn,
#                                            schema_trf=schema_trf,
#                                            otu_json=otu_jsonfi,
#                                            ingroup_mrca=ingroup_mrca)
#     filteredScrape = FilterBlast(data_obj, ids)
#     blast_dir = "tests/data/precooked/fixed/tte_blast_files"
#     filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
#     filteredScrape.remove_identical_seqs()
#     assert len(filteredScrape.new_seqs_otu_id) == 38

#     # fails with AttributeError: 'FilterBlast' object has no attribute 'mrca_ncbi'
