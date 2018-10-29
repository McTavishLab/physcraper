# ##############
# Martha Kandziora
# October 29th, 2017
# relabel tip names in phylogenetic trees
# define workdir and tree filename
# v. 0.1
# ##############


import pandas as pd
import sys


# ##################
# functions 

def read_otu_info(fn):
	df = pd.read_csv(fn, sep=',', header=None, index_col=False,
	                     names=[
	                         'otu_id',
	                         'tax_name',
	                         'gi',
	                         'accession',
	                         'last_blsated',
	                         'status',
	                         'ott_d',
	                         'ncbi_id',
	                         'GenBank info',
	                        
	                     ])
	return df



# ##################
# input

input_tr_fn = sys.argv[2]
input_workdir = sys.argv[1]


otu_fn = "{}/otu_seq_info.csv".format(input_workdir)
tr_fn = "{}/{}".format(input_workdir, input_tr_fn)
# ######
otu_info = read_otu_info(otu_fn)

with open(tr_fn, "r") as label_new:
	new_tree = label_new.read()
	print(new_tree)
	with open("{}_relabel".format(tr_fn), "wt") as fout:
		for index, row in otu_info.iterrows():
			print(row['otu_id'])
			ident = row['gi']
			if ident == "-":
				ident = row["otu_id"]

			new_tree = new_tree.replace("{}:".format(row['otu_id']), "{}_{}:".format(row['tax_name'].replace(" ", "_"), ident))

		fout.write(new_tree)
