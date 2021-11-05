import json
import physcraper
import sys
import urllib2
from Bio import Entrez




Entrez.email = 'ejmctavish@gmail.com'


configfi = "example.config"

conf = physcraper.ConfigObj(configfi)
workdir="ws-tests/tmp"

ids =  physcraper.IdDicts(conf, workdir=workdir)



otu_dict = json.loads(open("tests/data/tmp/otu_info.json").read())

#erase taxon info before search, for test
for item in otu_dict:
    try:
        del otu_dict[item]["^ncbi:taxon"]
        del otu_dict[item]["^ot:ottId"]
    except:
        pass

for item in otu_dict:
    acc = otu_dict[item].get("^ncbi:accession")
    if otu_dict[item].get("^ot:ottId"):
        sys.stdout.write("ott_id found\n")
        if not otu_dict[item].get("^ncbi:taxon"):
            otu_dict[item]["^ncbi:taxon"] = ids.ott_to_ncbi(otu_dict[item]["^ot:ottId"])
    elif acc:    
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, retmode="xml")
            read_handle = Entrez.read(handle)[0]
            tax_id = None
            for rec in read_handle['GBSeq_feature-table'][0]['GBFeature_quals']:
                if rec[u'GBQualifier_name'] == 'db_xref':
                    tax_id = rec[u'GBQualifier_value']
                    if tax_id.startswith('taxon'):
                        ncbi_id = tax_id.split(":")[1]
                if rec[u'GBQualifier_name'] == 'organism':
                    tax_name = u'GBQualifier_value'
            otu_dict[item]["^ncbi:taxon"] = ncbi_id
            try:
                otu_dict[item]["^ot:ottId"] = ids.ncbi_to_ott[int(ncbi_id)]
            except KeyError: 
                otu_dict[item]["^ot:ottId"] = None
                otu_dict[item]["ncbi:taxon"] = tax_name
        except:
            sys.stderr.write(err)
    else:
        sys.stderr.write("no taxon and no accession number\n")
