from physcraper import IdDicts, ConfigObj
import random

configfi = "tests/data/test.config"
workdir="tests/data/tmp/owndata"

def test_id_dicts():
    ids = IdDicts(configfi)
    selection  = random.sample(ids.ott_to_ncbi.keys(), 10)
    for ott_id in selection:
        ncbi_id = ids.ott_to_ncbi[ott_id]
        assert ids.ncbi_to_ott[ncbi_id] == ott_id
    assert(ids.config.hitlist_size == 100)

    

def test_id_no_config():
    ids = IdDicts()
    selection  = random.sample(ids.ott_to_ncbi.keys(), 10)
    for ott_id in selection:
        ncbi_id = ids.ott_to_ncbi[ott_id]
        assert ids.ncbi_to_ott[ncbi_id] == ott_id
    assert(ids.config.hitlist_size == 10)

test_id_dicts()