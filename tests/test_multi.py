import os
import json
import copy
import filecmp
from pytest import mark
from contextlib import contextmanager

from physcraper import OtuJsonDict, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape
from physcraper import opentree_helpers
#


seqaln= "tests/data/tiny_test_example/test.fas"
mattype="fasta"
trfn= "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"

workdir="tests/data/tmp/multi"

configfi = "tests/data/test.config"
id_to_spn = "tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)

"""Generates the files needed for the tests.
"""

if not os.path.exists("{}".format(workdir)):
    os.makedirs("{}".format(workdir))

conf_base = ConfigObj(configfi)

ids_base = IdDicts(configfi)

otu_json = OtuJsonDict(id_to_spn, ids_base)
with open(otu_jsonfi,"w") as outfile:
    json.dump(otu_json, outfile)

"""Extract ott ids from `otu_json` dictionary (???):
"""

ottids = [otu_json[ite]['^ot:ottId'] for ite in otu_json]


"""To get the MRCA of the matched taxa, use function `opentree_helpers.get_mrca_ott()`:
"""

mrca = opentree_helpers.get_mrca_ott(ottids)


"""Generate the alignment-taxon object:
"""

data_obj_base = generate_ATT_from_files(alnfile=seqaln,
                             aln_schema=mattype,
                             workdir=workdir,
                             configfile=configfi,
                             treefile=trfn,
                             tree_schema = schema_trf,
                             otu_json=otu_jsonfi,
                             search_taxon=mrca)


def test_otu_dict():
    assert(data_obj_base.otu_dict.keys()) == set(['2029_doronicum', 'S_doronicum', 'S_lagascanus', 'S_lopezii', 'S_scopolii'])




def test_add_all():
    threshold = 2
    conf = copy.deepcopy(conf_base)
    data_obj = copy.deepcopy(data_obj_base)
    ids = copy.deepcopy(ids_base)
    filteredScrape = PhyscraperScrape(data_obj, ids)
    filteredScrape._blasted = 1
    filteredScrape.config.spp_threshold = threshold
    filteredScrape.read_blast_wrapper(blast_dir="tests/data/precooked/fixed/tte_blast_files")
    filteredScrape.remove_identical_seqs()
    sp_d = filteredScrape.make_sp_dict(filteredScrape.new_seqs_otu_id)
    assert len(sp_d) == 4
    for taxon in sp_d:
        assert len(sp_d[taxon]) <= threshold



def test_no_mrca():
    search_taxon = None
    # setup the run
    if not os.path.exists("{}".format(workdir)):
        os.makedirs("{}".format(workdir))
    conf = copy.deepcopy(conf_base)
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)
    filteredScrape = PhyscraperScrape(data_obj, ids, search_taxon)
    filteredScrape.threshold = 5
    assert filteredScrape.mrca_ncbi == 795077
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    filteredScrape._blasted = 1
    filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
    filteredScrape.remove_identical_seqs()
    assert len(filteredScrape.new_seqs_otu_id) in [23,14] #Blurghhh, local vs remote searches get diffenrt number of seqs!



def test_remove_identical_seqs():
    conf = copy.deepcopy(conf_base)
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)
    # print("start")
    scraper = PhyscraperScrape(data_obj, ids)
    scraper._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    #scraper.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    scraper.read_blast_wrapper(blast_dir=blast_dir)
    #print scraper.ncbi_mrca
    assert(len(scraper.new_seqs) == 0)
    assert(len(scraper.data.aln) == 5)
    assert len(scraper.new_seqs_otu_id) == 14
    #Now that we are pulling the full remote sequences, we don'thave any identical seuqnces in the test.


def test_blocklist():
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)
    # print("start")
    scraper = PhyscraperScrape(data_obj, ids)
    scraper.blocklist.append('JX895264.1')
    scraper._blasted = 1
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    #scraper.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
    scraper.read_blast_wrapper(blast_dir=blast_dir)
    #print scraper.ncbi_mrca
    assert(len(scraper.new_seqs) == 0)
    assert(len(scraper.data.aln) == 5)
    assert len(scraper.new_seqs_otu_id) == 13
    aln_path1 = scraper.data.write_aln()
    scraper.align_new_seqs()
    assert len(scraper.data.aln) == 18
    os.remove(aln_path1)





localblast = mark.localblast

@localblast
def test_remove_taxa_aln_tre():
    conf = copy.deepcopy(conf_base)
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)

    filteredScrape =  PhyscraperScrape(data_obj, ids)

    len_aln_before = len(filteredScrape.data.aln.as_string('phylip'))
    len_tre_before = len(filteredScrape.data.tre.as_string(schema="newick"))
    namespace_before = len(filteredScrape.data.aln.taxon_namespace)
    namespace_tre_before = len(filteredScrape.data.tre.taxon_namespace)

    for tax in filteredScrape.data.aln.taxon_namespace:
        filteredScrape.data.remove_taxa_aln_tre(tax.label)
        break

    len_aln_after = len(filteredScrape.data.aln.as_string('phylip'))
    len_tre_after = len(filteredScrape.data.tre.as_string(schema="newick"))
    namespace_after = len(filteredScrape.data.aln.taxon_namespace)
    namespace_tre_after = len(filteredScrape.data.tre.taxon_namespace)

    assert len_aln_before != len_aln_after
    assert len_tre_before != len_tre_after
    assert namespace_before != namespace_after
    assert namespace_tre_before != namespace_tre_after


slow = mark.slow


#https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/


#@slow
def test_run_raxml():

    workdir = "tests/output/test_run_raxml"
    absworkdir = os.path.abspath(workdir)
    conf = copy.deepcopy(conf_base)
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)

    data_obj.workdir = absworkdir

    scraper = PhyscraperScrape(data_obj, ids)
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"
    scraper._blasted = 1

    # run needed functions
    # scraper.run_blast_wrapper()
    scraper.read_blast_wrapper(blast_dir=blast_dir)


    scraper.est_full_tree()
    scraper.calculate_final_tree(boot_reps=5)
    # scraper.generate_streamed_alignment()
    assert os.path.exists("{}/RAxML_bestTree.{}".format(scraper.rundir, scraper.date))
    assert os.path.exists("{}/updated_taxonname.tre".format(scraper.outputsdir))
    # scraper.generate_streamed_alignment()

    # if os.path.exists("{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date)):
    #   fn = "{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date)
    #   os.system("mv {}  ./previous_run/{}".format(fn, fn))

    # CWD = os.getcwd()
    # print(CWD)
    # os.chdir(scraper.workdir)
    # with timeout(30):
    #   print('entering block')
    #   scraper.calculate_bootstrap()

    #   # scraper.generate_streamed_alignment()
    #   import time
    #   time.sleep(50)
    #   print('This should never get printed because the line before timed out')
    # # os.chdir(CWD)
    # # scraper.calculate_boostrap()
    # assert os.path.exists("{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date))


def test_run_align():
    workdir = "tests/output/test_write_unaligned"
    absworkdir = os.path.abspath(workdir)
    conf = copy.deepcopy(conf_base)
    ids = copy.deepcopy(ids_base)
    data_obj = copy.deepcopy(data_obj_base)
    data_obj.workdir = absworkdir
    scraper = PhyscraperScrape(data_obj, ids)
    blast_dir = "tests/data/precooked/fixed/tte_blast_files"

    # run needed functions
    # scraper.run_blast_wrapper()
    scraper.read_blast_wrapper(blast_dir=blast_dir)

    assert(scraper.data.otu_dict['otuPS1']['^ncbi:taxon'] == int(scraper.data.otu_dict['2029_doronicum']['^ncbi:taxon']))
    scraper.data.aln.write(path="{}/myfilename.fas".format(scraper.workdir), schema='fasta')
    scraper.write_new_seqs(filename="only_new_seqs.aln")

    #scraper.est_full_tree()
    # scraper.generate_streamed_alignment()
    #assert os.path.exists("{}/RAxML_bestTree.{}".format(scraper.workdir, scraper.date))


def test_write_labelled():
    expected_tree_path = "tests/data/expected_output/labelled.tre"
    expected_aln_path = "tests/data/expected_output/labelled.fas"
    expected_tree_path_ottid = "tests/data/expected_output/labelled_ottid.tre"
    expected_aln_path_ottid = "tests/data/expected_output/labelled_ottid.fas"
    data_obj = copy.deepcopy(data_obj_base)
    treepath = 'tests/data/tmp/labelled_test.tre'
    alnpath = 'tests/data/tmp/labelled_test.fas'

    data_obj.write_labelled(label='^user:TaxonName', filename='labelled', direc='tests/data/tmp/', norepeats = False)

    assert os.path.isfile(treepath)
    assert os.path.isfile(alnpath)
    assert filecmp.cmp(treepath, expected_tree_path)
    assert filecmp.cmp(alnpath, expected_aln_path)
    os.remove(treepath)
    os.remove(alnpath)

    treepath_ottid = 'tests/data/tmp/labelled_ottid.tre'
    alnpath_ottid = 'tests/data/tmp/labelled_ottid.fas'

    data_obj.write_labelled(label='^ot:ottId', filename='labelled_ottid', direc='tests/data/tmp/', norepeats = False)

    assert os.path.isfile(treepath_ottid)
    assert os.path.isfile(alnpath_ottid)
    assert filecmp.cmp(treepath_ottid, expected_tree_path_ottid)
    assert filecmp.cmp(alnpath_ottid, expected_aln_path_ottid)
    os.remove(treepath_ottid)
    os.remove(alnpath_ottid)

    data_obj.write_otus(schema='table', direc='tests/data/tmp/')
    data_obj.write_otus(schema='json', direc='tests/data/tmp/')
    assert os.path.isfile('tests/data/tmp/otu_info_test.json')
    assert os.path.isfile('tests/data/tmp/otu_info_test.csv')
    os.remove('tests/data/tmp/otu_info_test.json')
    os.remove('tests/data/tmp/otu_info_test.csv')

