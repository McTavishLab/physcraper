from physcraper import IdDicts, ConfigObj
import random

configfi = "tests/data/test.config"
workdir="tests/data/tmp/owndata"

def test_id_dicts():
    conf = ConfigObj(configfi, interactive=True)
    ids = IdDicts(conf, workdir=workdir)
    selection  = random.sample(ids.ott_to_ncbi.keys(), 10)
    for ott_id in selection:
        ncbi_id = ids.ott_to_ncbi[ott_id]
        assert ids.ncbi_to_ott[ncbi_id] == ott_id
    


test_id_dicts()