import sys
from physcraper.concat import Concat

from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
concat = mark.concat



workdir_its = "tests/data/precooked/concat_pre"
workdir_ets = "tests/data/precooked/concat_pre"
email = "martha.kandziora@yahoo.com"
pickle_fn = "final_ATT_checkpoint.p"

workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": "its_{}".format(pickle_fn)}, 
            "ets": {"workdir": workdir_ets, "pickle": "ets_{}".format(pickle_fn)}}

# get to test status
@concat
def test():
    concat = Concat(workdir_comb, email)
    for item in genelist.keys():
        concat.load_single_genes(genelist[item]["workdir"], genelist[item]["pickle"], item)

    concat.combine()
    concat.sp_seq_counter()

    # print("tests if nothing gets lost from loading single runs to make sp_counter")

    counter = 0
    for sp in concat.sp_counter:
        for gene in concat.sp_counter[sp]:
            counter += concat.sp_counter[sp][gene]
    # print(counter)

    single_run_items = 0
    for item in concat.single_runs:
        single_run_items += (len(concat.single_runs[item].aln.taxon_namespace))

    assert counter == single_run_items
