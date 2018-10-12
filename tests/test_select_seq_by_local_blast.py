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
    conf = physcraper.ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = physcraper.IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()

filteredScrape =  physcraper.FilterBlast(data_obj, ids)
filteredScrape._blasted = 1
blast_dir = "tests/data/precooked/fixed/tte_blast_files"
# filteredScrape.gi_list_mrca = pickle.load(open("tests/data/precooked/gi_list_mrca.p", 'rb'))
filteredScrape.read_blast(blast_dir=blast_dir)
filteredScrape.remove_identical_seqs()
filteredScrape.sp_dict(downtorank)
filteredScrape.make_sp_seq_dict()

##this is the code of the first part of how many seq to keep. if threshold is bigger than number of seq for sp, just add all
# print("start test")


print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
count = 0
for tax_id in filteredScrape.sp_d:
    # print(filteredScrape.sp_d[tax_id])
    print(tax_id)
    # if len(filteredScrape.sp_d[tax_id]) > treshold:
    count_dict = filteredScrape.count_num_seq(tax_id)
    print(count_dict)
    if count_dict["new_taxon"]:
        print("new stuff")
        if count_dict["query_count"] < treshold:
            count += count_dict["query_count"]
        if count_dict["query_count"] > treshold:
            count += treshold
    if count_dict["new_taxon"] is False:
        print("old guy")
        if count_dict["query_count"] >= 1:
            if count_dict["seq_present"] < treshold:
                print("something can be added")
                count += treshold-count_dict["seq_present"]
            if count_dict["seq_present"] > treshold:
                count += 0
    print(count)
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

filteredScrape.how_many_sp_to_keep(treshold, selectby)

print("filteredScrape.filtered_seq")
print(filteredScrape.filtered_seq)
        # if tax_id in filteredScrape.sp_seq_d.keys():
        #     # print(tax_id)
        #     seq_present = count_dict["seq_present"]
        #     query_count = count_dict["query_count"]
        #     # for item in filteredScrape.sp_d[tax_id]:
        #     print(query_count, seq_present, treshold)
        #     if seq_present >= 1 and seq_present < treshold and count_dict["new_taxon"] == False and query_count != 0:
        #         if query_count + seq_present > treshold:
        #             print("loop for write blast files")
        #             taxonfn = filteredScrape.loop_for_write_blast_files(tax_id)
        #             for element in filteredScrape.sp_d[tax_id]:
        #                 if '^ot:ottTaxonName' in element:
        #                     blast_seq = "{}".format(element['^ot:ottTaxonName'])
        #                     blast_seq = blast_seq.replace(" ", "_")
        #                     blast_db = "{}".format(element['^ot:ottTaxonName'])
        #                     blast_db = blast_db.replace(" ", "_")
        #                     print(blast_db, blast_seq)
        #             if filteredScrape.downtorank is not None:
        #                 taxonfn = tax_id
        #             print(filteredScrape.data.workdir, taxonfn, taxonfn)
        #             local_blast.run_local_blast(filteredScrape.data.workdir, taxonfn, taxonfn)
        #             # print(filteredScrape.sp_seq_d[tax_id], taxonfn, treshold, seq_present)
        #             filteredScrape.select_seq_by_local_blast(filteredScrape.sp_seq_d[tax_id], taxonfn, treshold, seq_present)
        #     elif seq_present == 0 and count_dict["new_taxon"] == True and query_count>=1:
        #         print('elif')
        #         for item in filteredScrape.sp_d[tax_id]:
        #             if '^ncbi:accession' in item:
        #                 filteredScrape.data.add_otu(item['^ncbi:accession'], filteredScrape.ids)
        #         blast_seq = filteredScrape.sp_seq_d[tax_id].keys()[0]
        #         print(blast_seq)
        #         # print(some)
        #         if len(blast_seq.split(".")) >=2: # type(blast_seq) == int:
        #             str_db = str(tax_id)
        #         else:
        #             str_db = str(blast_seq)
        #         blast_db = filteredScrape.sp_seq_d[tax_id].keys()[1:]
        #         # write files for local blast first:
        #         seq = filteredScrape.sp_seq_d[tax_id][blast_seq]
        #         local_blast.write_blast_files(filteredScrape.data.workdir, str_db, seq) #blast qguy
        #         # print(blast_db)
        #         for blast_key in blast_db:
        #             seq = filteredScrape.sp_seq_d[tax_id][blast_key]
        #             local_blast.write_blast_files(filteredScrape.data.workdir, blast_key, seq, db=True, fn=str_db) #local db
        #         # make local blast of sequences
        #         if filteredScrape.downtorank is not None:
        #             str_db = tax_id
        #         local_blast.run_local_blast(filteredScrape.data.workdir, str_db, str_db)
        #         if len(filteredScrape.sp_seq_d[tax_id]) + seq_present >= treshold:
        #             filteredScrape.select_seq_by_local_blast(filteredScrape.sp_seq_d[tax_id], str_db, treshold, seq_present)
        #         elif len(filteredScrape.sp_seq_d[tax_id]) + seq_present < treshold:
        #             filteredScrape.add_all(tax_id)
        #     else:
        #         print("nothing will be added....to much is present")
        #     for otu in filteredScrape.data.otu_dict:
        #         # print(filteredScrape.data.otu_dict[otu])
        #         if filteredScrape.data.otu_dict[otu]['^ncbi:taxon'] == tax_id:
        #             print(filteredScrape.data.otu_dict[otu])

print(count, len(filteredScrape.filtered_seq) )
print(filteredScrape.filtered_seq.keys())
print(filteredScrape.sp_d)


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