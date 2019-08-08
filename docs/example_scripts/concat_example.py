
from physcraper import wrappers

workdir_ITS = "tests/data/PS_tiny_comb_its"  # folder of single-gene data to be concatenated
workdir_ETS = "tests/data/PS_tiny_comb_ets"  # folder of single-gene data to be concatenated
email = "mk@xy.zt"  # email address of user
percentage = 0.4  # percentage rof bases needed to be part of the concatenated data

num_threads = 4  # number of threads to use, to make it run faster

workdir_comb = "docs/example_scripts/output/nr_concat"  # folder to save the concatenated data
user_concat_fn = None  # if you want to define which sequences to be concatenated, provide file path here
backbone = False  # if you want to fix a backbone topology instead of recalculating whole relationship, set it to True

# dictionary that comprises the information of the different folder to be concatenated,
# add variables of your working directories and the name of the loci
genelist = {"ITS": {"workdir": workdir_ITS, "pickle": "final_ATT_checkpoint.p"},  
            "ETS": {"workdir": workdir_ETS, "pickle": "final_ATT_checkpoint.p"}
            }

conc = wrappers.concat(genelistdict=genelist, workdir_comb=workdir_comb,
                       email=email, num_threads=num_threads, percentage=percentage, user_concat_fn=user_concat_fn, backbone=backbone)