from __future__ import print_function, absolute_import 


expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'ncbi_parser_names_fn', 'ncbi_parser_nodes_fn', 'gb_id_filename', 'unmapped', 'hitlist_size', 'id_pickle', 'blast_loc', 'url_base', 'ott_ncbi', 'email', 'e_value_thresh', 'blastdb']


def test_config():
    from physcraper import ConfigObj
    configfi = "tests/data/test.config"
    conf = ConfigObj(configfi, interactive=False)
    assert conf.email == 'ejmctavish@gmail.com'
    assert conf.url_base == None
    assert conf.__dict__.keys() == expected_keys
