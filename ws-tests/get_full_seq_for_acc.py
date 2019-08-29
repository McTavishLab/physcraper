import physcraper


workdir = "tests/output/opentree"
configfi = "tests/data/remotencbi.config"




acc_ids = ['NG_058655.1', 'AF466085.1']
conf = physcraper.ConfigObj(configfi, interactive=False)
ids = physcraper.IdDicts(conf, workdir=workdir)

for gb_id in acc_ids:
    read_handle = self.ids.entrez_efetch(gb_id)
    seqlen = len(read_handle[0][u'GBSeq_sequence'])
