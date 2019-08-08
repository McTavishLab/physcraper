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
    configfi = "example.config"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    otu_jsonfi = "tests/data/tmp/owndata/otu_dict.json".format(workdir)


    conf = ConfigObj(configfi, interactive=False)

    data_obj = generate_ATT_from_files(seqaln=seqalnmiss, 
                                     mattype=mattype, 
                                    workdir=workdir,
                                     config_obj=conf,
                                     treefile=treefile,
                                     schema_trf = schema_trf,
                                     otu_json=otu_jsonfi,
                                     ingroup_mrca=None)

    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == '2029_doronicum':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"

    #----------------------------------------------------

    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                     mattype=mattype, 
                                    workdir=workdir,
                                     config_obj=conf,
                                     treefile=treefilemiss,
                                     schema_trf = schema_trf,
                                     otu_json=otu_jsonfi,
                                     ingroup_mrca=None)



    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == 'S_scopolii':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"




    #----------------------------------------------------


    aln = DnaCharacterMatrix.get(path=seqalnmiss, schema=mattype)

    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore UGH


    tre = Tree.get(path=treefile,
                   schema="newick",
                    preserve_underscores=True,
                    taxon_namespace=aln.taxon_namespace)



    assert aln.taxon_namespace == tre.taxon_namespace
    assert aln.taxon_namespace is tre.taxon_namespace


    treed_taxa = set()
    for leaf in tre.leaf_nodes():
        treed_taxa.add(leaf.taxon)
    aln_tax = set()
    for tax, seq in aln.items():
        aln_tax.add(tax)


    prune = treed_taxa ^ aln_tax

    assert len(prune) == 1
    assert list(prune)[0].label == '2029_doronicum'

    #----------------

    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

    assert aln.taxon_namespace
    for tax in aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")  # Forcing all spaces to underscore UGH


    tre = Tree.get(path=treefilemiss,
                   schema="newick",
                    preserve_underscores=True,
                    taxon_namespace=aln.taxon_namespace)



    assert aln.taxon_namespace == tre.taxon_namespace
    assert aln.taxon_namespace is tre.taxon_namespace


    treed_taxa = set()
    for leaf in tre.leaf_nodes():
        treed_taxa.add(leaf.taxon)
    aln_tax = set()
    for tax, seq in aln.items():
        aln_tax.add(tax)


    prune = treed_taxa ^ aln_tax

    assert len(prune) == 1
    assert list(prune)[0].label == 'S_scopolii'


    # ----------------------------



    seqaln= "tests/data/tiny_test_example/test.fas"
    seqalnmiss= "tests/data/tiny_test_example/test_missingseq.fas"
    mattype="fasta"
    treefile= "tests/data/tiny_test_example/test.tre"
    treefilemiss= "tests/data/tiny_test_example/test_missingtip.tre"
    schema_trf = "newick"
    workdir="tests/output/owndata"
    configfi = "example.config"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
    otu_jsonfi = "tests/data/tmp/owndata/otu_dict.json".format(workdir)

    data_obj = generate_ATT_from_files(seqaln=seqalnmiss, 
                                     mattype=mattype, 
                                     workdir=workdir,
                                     config_obj=conf,
                                     treefile=treefilemiss,
                                     schema_trf = schema_trf,
                                     otu_json=otu_jsonfi,
                                     ingroup_mrca=None)




    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == '2029_doronicum':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"

    for otu in data_obj.otu_dict:
        if data_obj.otu_dict[otu][u'^ot:originalLabel'] == 'S_scopolii':
            assert data_obj.otu_dict[otu]['^physcraper:status'] == "deleted in reconciliation"

