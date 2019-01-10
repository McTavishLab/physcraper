import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, ConfigObj, IdDicts
#

# Setup runs for the concat test functions

# tiny its
seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"
workdir = "tests/data/PS_tiny_comb_its"
configfi = "tests/data/localblast.config"
threshold = 2
selectby = "blast"



wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         id_to_spn,
                         configfi,
                         selectby=selectby)

# tiny ets
seqaln = "tests/data/tiny_comb_ets/tiny_comb_ets.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_ets/tiny_comb_ets.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_comb_ets/nicespl.csv"

workdir = "tests/data/PS_tiny_comb_ets"
configfi = "tests/data/localblast.config"
treshold = 2
selectby = "blast"



wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         id_to_spn,
                         configfi,
                         selectby=selectby)
