import os
import json
from dendropy import Tree, \
      DnaCharacterMatrix, \
      DataSet, \
      datamodel
from physcraper import ConfigObj, generate_ATT_from_files, AlignTreeTax, OtuJsonDict


def test_reconcile():
    #------------------------
    seqaln= "tests/data/tiny_test_example/test.fas"
    seqalnmiss= "tests/data/tiny_test_example/test_missingseq.fas"
    mattype="fasta"
    treefile= "tests/data/tiny_test_example/test.tre"
    treefilemiss= "tests/data/tiny_test_example/test_missingtip.tre"
    schema_trf = "newick"
    workdir="tests/output/owndata"
    configfi = "docs/examples/example.config"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    otu_jsonfi = "tests/data/tiny_test_example/otu_dict.json"


    data_obj = generate_ATT_from_files(alnfile=seqalnmiss, 
                                     aln_schema=mattype, 
                                     workdir=workdir,
                                     configfile=configfi,
                                     treefile=treefile,
                                     tree_schema=schema_trf,
                                     otu_json=otu_jsonfi,
                                     search_taxon=None)

    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == '2029_doronicum':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"

    #----------------------------------------------------

    data_obj = generate_ATT_from_files(alnfile=seqaln, 
                                     aln_schema=mattype, 
                                     workdir=workdir,
                                     configfile=configfi,
                                     treefile=treefilemiss,
                                     tree_schema=schema_trf,
                                     otu_json=otu_jsonfi,
                                     search_taxon=None)


    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == 'S_scopolii':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"


    # ----------------------------



    seqalnmiss= "tests/data/tiny_test_example/test_missingseq.fas"
    treefilemiss= "tests/data/tiny_test_example/test_missingtip.tre"

    data_obj = generate_ATT_from_files(alnfile=seqalnmiss, 
                                       aln_schema=mattype, 
                                       workdir=workdir,
                                       configfile=configfi,
                                       treefile=treefilemiss,
                                       tree_schema = schema_trf,
                                       otu_json=otu_jsonfi,
                                      search_taxon=None)




    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == '2029_doronicum':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == 'S_scopolii':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"

