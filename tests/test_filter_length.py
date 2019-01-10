import os
import json
import sys
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts, generate_ATT_from_files, FilterBlast
#
from pytest import mark

localblast = mark.localblast

@localblast
def test_filter_length():
    seqaln = "tests/data/tiny_test_example/test.fas"
    mattype = "fasta"
    trfn = "tests/data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"

    workdir = "tests/output/test_selectbylength"
    configfi = "tests/data/test.config"
    otu_jsonfi = "{}/otu_dict.json".format(workdir)
    threshold = 2
    selectby = "length"
    downtorank = "species"
    add_unpubl_seq = None
    blacklist=None
                    
    id_to_spn_addseq_json=None
    ingroup_mrca=None
    shared_blast_folder=None


    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))

    conf = ConfigObj(configfi)
    ids = IdDicts(conf, workdir=workdir)

    
    otu_json = OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi, "w"))

      
#            sync_names()
    sys.stdout.write("setting up Data Object\n")
    sys.stdout.flush()
    #read the config file into a configuration object
    conf = ConfigObj(configfi)

    #Generate an linked Alignment-Tree-Taxa object
    data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                           mattype=mattype,
                                           workdir=workdir,
                                           config_obj=conf,
                                           treefile=trfn,
                                           schema_trf=schema_trf,
                                           otu_json=otu_jsonfi,
                                           ingroup_mrca=ingroup_mrca)

    # Prune sequnces below a certain length threshold
    # This is particularly important when using loci that have been de-concatenated,
    # as some are 0 length which causes problems.
    data_obj.prune_short()
    data_obj.write_files()
    data_obj.write_labelled(label="^ot:ottTaxonName", add_gb_id=True)
    data_obj.write_otus("otu_info", schema="table")
    data_obj.dump()
    sys.stdout.write("setting up id dictionaries\n")
    sys.stdout.flush()
    ids = IdDicts(conf, workdir=workdir, mrca=ingroup_mrca)

    # Now combine the data, the ids, and the configuration into a single physcraper scrape object
    filteredScrape = FilterBlast(data_obj, ids)
    filteredScrape.add_setting_to_self(downtorank, threshold)
    filteredScrape.blacklist = blacklist
   
    sys.stdout.write("BLASTing input sequences\n")
    if shared_blast_folder:
        filteredScrape.blast_subdir = shared_blast_folder
    else:
        shared_blast_folder = None
    filteredScrape.run_blast_wrapper()
    filteredScrape.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
    filteredScrape.remove_identical_seqs()
    filteredScrape.dump()
    sys.stdout.write("Filter the sequences\n")
    length_unfiltered = len(filteredScrape.new_seqs)

    if threshold is not None:
        filteredScrape.sp_dict(downtorank)
        filteredScrape.make_sp_seq_dict()
        filteredScrape.how_many_sp_to_keep(selectby=selectby)
        filteredScrape.replace_new_seq()

    length_filtered = len(filteredScrape.new_seqs)

    assert length_filtered != length_unfiltered