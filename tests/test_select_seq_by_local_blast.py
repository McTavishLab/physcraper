import sys
import os
#from physcraper import ConfigObj, IdDicts, FilterBlast
import pickle#
import physcraper
import physcraper.local_blast as local_blast


sys.stdout.write("\ntests select_seq_by_local_blast\n")

# tests select_seq_local_blast_test
# tests if the building of select_seq_from local blast is selecting the right amount of species
workdir = "tests/output/test_select_seq_local_blast"
configfi = "tests/data/test.config"
treshold = 2
selectby = "blast"
downtorank = None

absworkdir = os.path.abspath(workdir)


try:
    conf = physcraper.ConfigObj(configfi, interactive=False)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  physcraper.FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
# filteredScrape.acc_list_mrca = pickle.load(open("tests/data/precooked/acc_list_mrca.p", 'rb'))
filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict()

##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
# print("start test")
count = 0
for tax_id in filteredScrape.sp_d:
    count_dict = filteredScrape.count_num_seq(tax_id)
    if count_dict["new_taxon"]:
        if count_dict["query_count"] < treshold:
            count += count_dict["query_count"]
        if count_dict["query_count"] > treshold:
            count += treshold
    if count_dict["new_taxon"] is False:
        if count_dict["query_count"] >= 1:
            if count_dict["seq_present"] < treshold:
                count += treshold-count_dict["seq_present"]
            if count_dict["seq_present"] > treshold:
                count += 0
filteredScrape.how_many_sp_to_keep(treshold, selectby)

try:
    assert count == len(filteredScrape.filtered_seq) and count>0
    sys.stdout.write("\ntest passed\n")
except:
    sys.stderr.write("\ntest failed\n")

# #added before
# #[429489224, 429489233, 429489188]
# {'^ncbi:taxon': 1268591, '^ncbi:title': 'Senecio scopolii subsp. scopolii clone JC4715-6 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence', '^ot:ottTaxonName': 'Senecio_scopolii_subsp._scopolii', '^physcraper:status': 'query', '^ot:ottId': 114544, '^ncbi:accession': 'JX895389.1', '^ncbi:gi': 429489224, '^physcraper:last_blasted': '1800/01/01'}
# {'^ncbi:taxon': 1268591, '^ncbi:title': 'Senecio scopolii subsp. scopolii clone JC4715-15 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence', '^ot:ottTaxonName': 'Senecio_scopolii_subsp._scopolii', '^physcraper:status': 'query', '^ot:ottId': 114544, '^ncbi:accession': 'JX895398.1', '^ncbi:gi': 429489233, '^physcraper:last_blasted': '1800/01/01'}
# {'^ncbi:taxon': 1268580, '^ncbi:title': 'Senecio lagascanus clone JC5600-6 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence', '^ot:ottTaxonName': 'Senecio_lagascanus', '^physcraper:status': 'query', '^ot:ottId': 640718, '^ncbi:accession': 'JX895353.1', '^ncbi:gi': 429489188, '^physcraper:last_blasted': '1800/01/01'}


# [u'JX895398.1', u'JX895353.1', u'JX895392.1', 'JX895513.1', 'JX895264.1']

# ## now only one scopolii
# 1268590: [{'^ncbi:taxon': 1268590, '^ncbi:title': 'Senecio scopolii subsp. floccosus 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence', '^ot:ottTaxonName': 'Senecio_scopolii_subsp._floccosus', '^physcraper:status': 'query', '^ot:ottId': 114541, '^ncbi:accession': 'JX895513.1', '^ncbi:gi': 429489348, '^physcraper:last_blasted': '1800/01/01'}],
# 1268581: {'^ncbi:taxon': 1268581, '^ncbi:title': 'Senecio lopezii clone JC3604-12 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence', '^ot:ottTaxonName': 'Senecio_lopezii', '^physcraper:status': 'query', '^ot:ottId': 688688, '^ncbi:accession': 'JX895264.1', '^ncbi:gi': 429489099, '^physcraper:last_blasted': '1800/01/01'}