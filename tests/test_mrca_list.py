# package import
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts

# define here your files
def test_mrca_list():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    workdir = "tests/output/impls_mrcalist_local"
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
