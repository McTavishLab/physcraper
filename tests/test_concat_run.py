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

	workdir_ITS = "tests/data/PS_tiny_comb_its"
	workdir_ETS = "tests/data/PS_tiny_comb_ets"
	email = "mk@xy.zt"
	percentage = 0.4

	pickle_fn = "final_ATT_checkpoint.p"

	workdir_comb = "./tests/output/concat_nr"
	genelist = {"ITS": {"workdir": workdir_ITS, "pickle": pickle_fn}, 
	            "ETS": {"workdir": workdir_ETS, "pickle": pickle_fn}
	            }

	conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
	                       email=email, percentage=percentage, user_concat_fn=None)

	assert isinstance(conc, Concat)
	print(os.path.join(workdir_comb, "physcraper"))
	shutil.rmtree(os.path.join(workdir_comb, "physcraper_runcopy"))