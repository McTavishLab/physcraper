from physcraper.concat import Concat

#
workdir_its = "data/tests/tiny_comb_its"
workdir_ets = "data/tests/tiny_comb_ets"
email = "martha.kandziora@yahoo.com"

pickle_fn = "scrape_checkpoint.p"

workdir_comb = "tests/output/impl_concat"
genelist = {"its": {"workdir": workdir_its, "pickle": pickle_fn}, "ets": {"workdir": workdir_ets, "pickle": pickle_fn}}

##############
# print("{}/{}".format(workdir, pickle_fn))
# scrape = pickle.load(open("{}/{}".format(workdir, pickle_fn),'rb'))
# print(scrape)
#############
# get to test status

concat = Concat(workdir_comb, email)
for item in genelist.keys():
    concat.load_single_genes(genelist[item]['workdir'], genelist[item]["pickle"], item)

concat.combine()
concat.sp_seq_counter()
sp_to_keep = concat.sp_to_keep()

print(sp_to_keep.keys())

print("tests sp_to_keep")

counter = 0
sp_keep = []
for sp in concat.sp_counter:
    for gene in concat.sp_counter[sp]:
        if concat.sp_counter[sp][gene] == 0:
            sp_keep.append(sp)

print(sp_keep)
print(sp_to_keep.keys())
print(len(sp_keep))

print(len(sp_to_keep.keys()))

try:
    assert set(sp_to_keep.keys()) == set(sp_keep)
    print("tests passed")
except:
    print("test fails")

