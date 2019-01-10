from __future__ import print_function, absolute_import 


#expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'ncbi_parser_names_fn', 'ncbi_parser_nodes_fn', 'gb_id_filename', 'unmapped', 'hitlist_size', 'id_pickle', 'blast_loc', 'url_base', 'ott_ncbi', 'email', 'e_value_thresh', 'blastdb']
# expected_keys = ['email', 'e_value_thresh', 'hitlist_size', 'blast_loc', 'blastdb', 'url_base',  'ncbi_parser_nodes_fn', 'ncbi_parser_names_fn', 'num_threads', 'gb_id_filename', 'delay',  'unmapped', 'seq_len_perc', 'seq_len_perc', 'trim_perc', 'phylesystem_loc', 'ott_ncbi', 'id_pickle' ]



def test_config():
    from physcraper import ConfigObj
    configfi = "tests/data/test.config"
    conf = ConfigObj(configfi, interactive=False)
    print(conf.__dict__.keys())


    if conf.blast_loc != "remote":
    	expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'ncbi_parser_names_fn', 'ncbi_parser_nodes_fn',  'maxlen', 'hitlist_size', 'gb_id_filename', 'delay', 'unmapped', 'trim_perc', 'url_base', 'ott_ncbi', 'blast_loc', 'id_pickle', 'email', 'e_value_thresh', 'blastdb']
    else:
		expected_keys = ['seq_len_perc', 'num_threads', 'phylesystem_loc', 'maxlen', 'hitlist_size', 'gb_id_filename', 'delay', 'unmapped', 'trim_perc', 'url_base', 'ott_ncbi', 'blast_loc', 'id_pickle', 'email', 'e_value_thresh']

    assert len(conf.email.split("@")) == 2
    assert conf.url_base == None
    assert conf.__dict__.keys() == expected_keys
