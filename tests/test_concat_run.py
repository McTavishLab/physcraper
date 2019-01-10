import os
import shutil
from physcraper import wrappers
from physcraper.concat import Concat

from pytest import mark
# you can actually do whatever
# ruftrum = mark.ruftrum will work and create a "ruftrum" test. 
# concatfull = mark.concatfull



# @concatfull
@mark.xfail
def test_concat_run():
    workdir_its = "tests/data/precooked/concat_pre"
    workdir_ets = "tests/data/precooked/concat_pre"
    email = "martha.kandziora@yahoo.com"
    pickle_fn = "final_ATT_checkpoint.p"
    email = "mk@xy.zt"
    percentage = 0.4
    workdir_comb = "./tests/output/concat_nr"
    genelist = {"ITS": {"workdir": workdir_its, "pickle": pickle_fn}, 
	            "ETS": {"workdir": workdir_ets, "pickle": pickle_fn}
	            }
    conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
	                       email=email, percentage=percentage, user_concat_fn=None)
    assert isinstance(conc, Concat)    
    print(os.path.join(workdir_comb, "physcraper"))    
    shutil.rmtree(os.path.join(workdir_comb, "physcraper_runcopy"))


def test_concat_combine():   

	workdir_its = "tests/data/precooked/concat_pre"
	workdir_ets = "tests/data/precooked/concat_pre"
	email = "martha.kandziora@yahoo.com"
	pickle_fn = "final_ATT_checkpoint.p"
	email = "mk@xy.zt"
	percentage = 0.4
	num_threads = 2
	workdir_comb = "./tests/output/concat_nr"
	genelistdict = {"its": {"workdir": workdir_its, "pickle": "its_{}".format(pickle_fn)}, 
	        "ets": {"workdir": workdir_ets, "pickle": "ets_{}".format(pickle_fn)}}
	conc = Concat(workdir_comb, email)
	for item in genelistdict.keys():
	    conc.load_single_genes(genelistdict[item]["workdir"], genelistdict[item]["pickle"], item)
	conc.combine()
	assert os.path.exists("{}/load_single_data.p".format(workdir_comb))
	   

	conc.sp_seq_counter()
	conc.get_largest_tre()
	assert conc.tre_start_gene == "ets"
	conc.make_sp_gene_dict()
	conc.make_alns_dict()

	for gene in conc.aln_all:
		conc.li((conc.aln_all[gene].taxon_namespace))
		assert len(conc.aln_all[gene]) == 20

	conc.concatenate_alns()
	conc.get_short_seq_from_concat(percentage)
	conc.remove_short_seq()
	assert len(conc.concatenated_aln) == 5

