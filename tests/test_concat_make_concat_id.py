import sys
import os
from physcraper.concat import Concat

#
workdir_its = "runs/tiny_comb_its"
workdir_ets = "runs/tiny_comb_ets"
email = "martha.kandziora@yahoo.com"
pickle_fn = "scrape_checkpoint.p"

workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, 
            "ets": {"workdir": workdir_ets, "pickle": pickle_fn}}


sys.stdout.write("\ntests Concat func make_concat_id_dict\n")


concat = Concat(workdir_comb, email)
print(os.getcwd())
for item in genelist.keys():
    concat.load_single_genes(genelist[item]['workdir'], genelist[item]["pickle"], item)

concat.combine()
spnl = []
for genename in concat.single_runs:
    for otu in concat.single_runs[genename].data.otu_dict.keys():
        data = concat.single_runs[genename].data.otu_dict[otu]
       
        if '^ot:ottTaxonName' in data:
            spn = concat.get_taxon_info('^ot:ottTaxonName', data)
        elif '^user:TaxonName' in data:
            spn = concat.get_taxon_info('^user:TaxonName', data)
        spnl.append(spn)

len_single = len(set(spnl))
# print(set(spnl))
# print(concat.sp_gi_comb)
len_concat_id_dict = len(concat.sp_gi_comb)


# print(len_single, len_concat_id_dict)


try:
    assert len_single == len_concat_id_dict
    print("tests passed")
except:
    print("test fails")
