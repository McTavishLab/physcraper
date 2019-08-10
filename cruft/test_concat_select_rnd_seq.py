
from copy import deepcopy
from physcraper.concat import Concat
from physcraper.helpers import debug
import sys

from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
concat = mark.concat


@concat
def test():
    workdir_its = "tests/data/precooked/concat_pre"
    workdir_ets = "tests/data/precooked/concat_pre"
    email = "martha.kandziora@yahoo.com"
    pickle_fn = "final_ATT_checkpoint.p"

    workdir_comb = "tests/output/impl_concat"
    genelist = {"its": {"workdir": workdir_its, "pickle": "its_{}".format(pickle_fn)}, 
            "ets": {"workdir": workdir_ets, "pickle": "ets_{}".format(pickle_fn)}}

    # get to test status

    sys.stdout.write("\ntests Concat func select_rnd_seq\n")


    concat = Concat(workdir_comb, email)
    for item in genelist.keys():
        concat.load_single_genes(genelist[item]['workdir'], genelist[item]["pickle"], item)


    concat.combine()
    concat.sp_seq_counter()
    sp_to_keep = concat.sp_to_keep()
    concat.get_largest_tre()


    # print("test: select rnd seq")
    count = 2
    concat.tmp_dict = deepcopy(concat.sp_acc_comb)
    # print("while")
    # part of make_sp_gene_dict

    len_before = len(concat.comb_seq)

    while len(concat.tmp_dict.keys()) >= 1:
        del_gi = {}
        for spn in concat.tmp_dict.keys():
            # print(spn)
            sp_to_keep_list = sp_to_keep.keys()
            # debug(sp_to_keep_list)
            if spn.replace(" ", "_") in sp_to_keep_list:
                tmp_gene = deepcopy(concat.genes_present)
                for gene in concat.tmp_dict[spn]:

                    tmp_gene.remove(gene)
                    # print("select_rnd_seq")
                    # print(spn, gene,del_gi, count)
                    del_gi = concat.select_rnd_seq(spn, gene, del_gi)
                    # print("now it should break")
                    break
                break
            break
        break

    len_after = len(concat.comb_seq[gene])

    #

    assert len_before + 1 == len_after
