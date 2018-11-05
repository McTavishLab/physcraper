from __future__ import print_function, absolute_import 


expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'blast_loc', 'gifilename', 'get_ncbi_taxonomy', 'hitlist_size', 'id_pickle', 'ott_ncbi', 'url_base', 'ncbi_dmp', 'email', 'e_value_thresh', 'blastdb']

def test_config():
    from physcraper import ConfigObj
    configfi = "tests/data/localblast.config"
    conf = ConfigObj(configfi, interactive=False)
    assert conf.email == 'ejmctavish@gmail.com'
    assert conf.url_base == None
    assert conf.__dict__.keys() == expected_keys
