from physcraper import wrappers


workdir_its = "./tiny_comb_its"
workdir_ets = "./tiny_comb_ets"
email = "martha.kandziora@yahoo.com"
percentage = 0.4  # how much missing data between datasets is allowed. Needs to be low, if you want to allow to incl sequences that are only found in the shorter of the two alignments.

pickle_fn = "scrape_checkpoint.p"


workdir_comb = "./tests/output/impl_concat_ownfile"
genelistdict = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, 
                "ets": {"workdir": workdir_ets, "pickle": pickle_fn},
                }

# Following file determines which sequences are concatenated together
user_concat_fn = "concatenation_input.csv"  # file must be in wd!

conc = wrappers.concat(genelistdict=genelistdict, workdir_comb=workdir_comb, 
					   email=email, percentage=percentage, user_concat_fn=user_concat_fn)
