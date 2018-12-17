import pickle
import sys
import os
import json
import physcraper 

sys.stdout.write("\ntests prune_short\n")

seqaln= "tests/data/tiny_test_example/test_extrashortseq.fas"
mattype="fasta"
treefile= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir="tests/output/test_trim"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)


def test_prune_short():
    if not os.path.exists("{}".format(workdir)):
            os.makedirs("{}".format(workdir))

    conf = physcraper.ConfigObj(configfi, interactive=False)
    ids = physcraper.IdDicts(conf, workdir=workdir)

    if os.path.exists(otu_jsonfi):
        otu_json = json.load(open(otu_jsonfi))
    else:
        otu_json = physcraper.OtuJsonDict(id_to_spn, ids)
        json.dump(otu_json, open(otu_jsonfi,"w"))


    data_obj = physcraper.generate_ATT_from_files(seqaln=seqaln, 
                                     mattype=mattype, 
                                     workdir=workdir,
                                     config_obj=conf,
                                     treefile=treefile,
                                     schema_trf = schema_trf,
                                     otu_json=otu_jsonfi,
                                     ingroup_mrca=None)


    data_obj.config.seq_len_perc = 0.9
    len_before = len(data_obj.tre.taxon_namespace)
    data_obj.prune_short()
    len_after = len(data_obj.tre.taxon_namespace)
    assert len_before > len_after
