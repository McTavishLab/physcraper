
import sys
from physcraper import ConfigObj

expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'blast_loc', 'get_ncbi_taxonomy', 'hitlist_size', 'id_pickle', 'ott_ncbi', 'url_base', 'ncbi_dmp', 'email', 'e_value_thresh', 'blastdb']

sys.stdout.write("\nTesting configuration object contents\n")
try:
    configfi = "tests/data/localblast.config"
    conf = ConfigObj(configfi)
    assert conf.email == 'ejmctavish@gmail.com'
    assert conf.url_base == None
    assert conf.__dict__.keys() == expected_keys
    sys.stdout.write("\nTest configread.py passed\n")
except:
	sys.stdout.write("\nTest configread.py FAILED (expected, test in progress)\n")