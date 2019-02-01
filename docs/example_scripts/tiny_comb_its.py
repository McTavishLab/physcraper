import sys
import os
import json
from physcraper import wrappers, OtuJsonDict, IdDicts, ConfigObj

# tiny ets
seqaln = "tests/data/tiny_comb_its/tiny_comb_its.fasta"
mattype = "fasta"
trfn = "tests/data/tiny_comb_its/tiny_comb_its.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/tiny_comb_its/nicespl.csv"

workdir = "docs/example_scripts/output/tiny_comb_its"
configfi = "tests/data/localblast.config"
threshold = 2
selectby = "blast"
downtorank = None


wrappers.filter_data_run(seqaln,
                         mattype,
                         trfn,
                         schema_trf,
                         workdir,
                         threshold,
                         id_to_spn,
                         configfi,
                         selectby=selectby,
                         downtorank=downtorank)
