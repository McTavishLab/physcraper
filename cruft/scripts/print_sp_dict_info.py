import sys
import os
import pickle  #
import physcraper 
import csv


scrape = pickle.load(open("/home/blubb/sync-TP-T470s/physcraper_MS/cluster_run1029/output/species/Senecioneae_ets_expand/scrape_checkpoint.p", 'rb'))


sp_d = scrape.sp_dict('species')
sp_info = {}


for k in sp_d:
    sp_info[k] = len(sp_d[k])
print(len(sp_d))
with open('taxon_sampling_spd.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in sp_info.items():
    	spn = scrape.ids.ncbi_parser.get_name_from_id(key)
        writer.writerow([key, spn, value])


sp_d = scrape.sp_dict('species')
sp_info = {}
print(len(sp_d))

for k in sp_d:
    sp_info[k] = len(sp_d[k])

with open('taxon_sampling_spd.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in sp_info.items():
    	spn = scrape.ids.ncbi_parser.get_name_from_id(key)
        writer.writerow([key, spn, value])
