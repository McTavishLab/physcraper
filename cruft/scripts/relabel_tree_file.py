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

	                         "otuID",
							"^ncbi:gi",
							"^ncbi:accession",
							"^ot:originalLabel",
							"^physcraper:last_blasted",
							"^physcraper:status",
							"^physcraper:TaxonName",
							"^ncbi:title",
							"^ncbi:taxon",
							"^ncbi:TaxonName",
							"^ot:ottId",
							"^ot:ottTaxonName"
                        
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
			print(row['otuID'])
			ident = row['^ncbi:accession']
			if ident == "-":
				ident = row["otuID"]

			new_tree = new_tree.replace("{}:".format(row['otuID']), "{}_{}:".format(row['^physcraper:TaxonName'].replace(" ", "_"), ident))

		fout.write(new_tree)
