import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, generate_ATT_from_files, AlignTreeTax, ConfigObj, IdDicts, FilterBlast
import pickle#

# tests if identical sequences are removed.
#
seqaln = "tests/data/tiny_test_example/test.fas"
mattype = "fasta"
trfn = "tests/data/tiny_test_example/test.tre"
schema_trf = "newick"
workdir = "tests/output/test_remove_id_seq"
configfi = "tests/data/test.config"
id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
add_local_seq = None
id_to_spn_addseq_json = None
blast_dir = "tests/data/precooked/fixed/tte_blast_files"


absworkdir = os.path.abspath(workdir)




# conf = ConfigObj(configfi)
# data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))

# for tax in data_obj.aln.taxon_namespace:
#     print(data_obj.aln[tax].symbols_as_string())

# # print(print(.aln.symbols_as_string())
# print(some)
try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    # print(data_obj.aln.symbols_as_string())
    # print(some)
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb" ))
except:
    # sys.stderr.write("run 'python tests/testfilesetup.py' to setup data files for tests. EXITING")
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
# filteredScrape.read_blast(blast_dir=blast_dir)

#############################

id_seq = ["TCGAAACCTGCATAGCAGAACGACCT-GTGAACATGTAAAAACAATTGGG-TGTTCTAAGTATCGGGCTCTTGTTCGATTTCTA-GGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGT-CTAAGGACGTCACGTCGACG-CAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGC--TT-GTTCCATGCATT--GCCGTT--CGCGGTGATTGCATTGAAACTTGCTTCTTTATAA-TTCATAAACGACTCTCGG-CAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCC-GAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGCGTCGCCC-CCATCAC---ACCTCTT-GACGGGGATGTTTGAATGGGGA-CGGAGATTGGTCTCCCGTTCCT---AAGGTGCGGTTGCCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCT--------------TATCGAGTTGTGTG--TTCCAAGAAGTAA-GGAATATCTCTTTAACGACCC-TAAAGTGTTGTCTCATG-ACGATGCTTCGACTGC",
            "TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC",
            "TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCGCGCGCGC",
            "TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCGCGCGCGC"
            ]
#############################
print("start test")

tmp_dict = dict((taxon.label, filteredScrape.data.aln[taxon].symbols_as_string()) for taxon in filteredScrape.data.aln)
print(tmp_dict)
old_seqs = tmp_dict.keys()
#Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
avg_seqlen = sum(filteredScrape.data.orig_seqlen)/len(filteredScrape.data.orig_seqlen) #HMMMMMMMM
assert filteredScrape.config.seq_len_perc <= 1
seq_len_cutoff = avg_seqlen*filteredScrape.config.seq_len_perc
# for gi, seq in filteredScrape.new_seqs.items():
count=1

for item in id_seq:
    print(item)
    gi =  1061375300
    if len(item.replace("-", "").replace("N", "")) > seq_len_cutoff:
        # otu_id = filteredScrape.data.add_otu(gi, filteredScrape.ids)
        # print(item)
      
        ott = "OTT_{}".format(count)
        count += 1

        otu_id = ott
        filteredScrape.data.otu_dict[otu_id] = {}
        filteredScrape.data.otu_dict[otu_id]['^ncbi:gi'] = gi
        filteredScrape.data.otu_dict[otu_id]['^ncbi:accession'] =   "KX494441"
        filteredScrape.data.otu_dict[otu_id]['^ncbi:title'] = "some random title"
        filteredScrape.data.otu_dict[otu_id]['^ncbi:taxon'] = "ncbi_id"
        filteredScrape.data.otu_dict[otu_id]['^ot:ottId'] = ott
        filteredScrape.data.otu_dict[otu_id]['^physcraper:status'] = "query"
        filteredScrape.data.otu_dict[otu_id]['^ot:ottTaxonName'] = "spn"
        filteredScrape.data.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"

        # otu_id = "Senecio_doronicum"
        filteredScrape.seq_dict_build(item, otu_id, tmp_dict)
for tax in old_seqs:
    try:
        del tmp_dict[tax]
    except KeyError:
        pass

print("end of seq dict build")
print(tmp_dict)
filteredScrape.new_seqs_otu_id = tmp_dict


# if len(filteredScrape.new_seqs) > 0:
#     filteredScrape.remove_identical_seqs()
#     filteredScrape.data.write_files() #should happen before aligning in case of pruning
#     if len(filteredScrape.new_seqs_otu_id) > 0:#TODO rename to something more intutitive
#         filteredScrape.write_query_seqs()
#         filteredScrape.align_query_seqs()
#         filteredScrape.data.reconcile()
#         filteredScrape.place_query_seqs()
#         filteredScrape.est_full_tree()

##################3
expected_add =1
print(len(filteredScrape.new_seqs_otu_id))

try:
    assert expected_add == len(filteredScrape.new_seqs_otu_id)
    print("test passed")
    print("todo: add check that newly added seq are checked. they are, but there is no test")
except:
    print("test failed")