import sys
from physcraper.concat import Concat

#
workdir_its = "tests/data/precooked/concat_pre"
workdir_ets = "tests/data/precooked/concat_pre"
email = "martha.kandziora@yahoo.com"
pickle_fn = "final_ATT_checkpoint.p"

workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": "its_{}".format(pickle_fn)}, 
            "ets": {"workdir": workdir_ets, "pickle": "ets_{}".format(pickle_fn)}}


from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
concat = mark.concat


@concat
def test():

    concat = Concat(workdir_comb, email)
    for item in genelist.keys():
        concat.load_single_genes(genelist[item]['workdir'], genelist[item]["pickle"], item)

    concat.combine()
    spnl = []
    for genename in concat.single_runs:
        for otu in concat.single_runs[genename].otu_dict.keys():
            data = concat.single_runs[genename].otu_dict[otu]
           
            if '^ot:ottTaxonName' in data:
                spn = concat.get_taxon_info('^ot:ottTaxonName', data)
            elif '^user:TaxonName' in data:
                spn = concat.get_taxon_info('^user:TaxonName', data)
            spnl.append(spn)

    len_single = len(set(spnl))
    # print(set(spnl))
    # print(concat.sp_acc_comb)
    len_concat_id_dict = len(concat.sp_acc_comb)


    assert len_single == len_concat_id_dict
