import pytest

def pytest_addoption(parser):

    parser.addoption("--runslow", action="store_true",
        help="run slow tests")
    parser.addoption("--runconcat", action="store_true",
        help="run concat tests")
    # parser.addoption("--runconcatfull", action="store_true",
    #     help="run concat full tests")
    parser.addoption("--runweb", action="store_true",
        help="run web tests")
    parser.addoption("--nolocalblast", action="store_true",
        help="do not run tests that rely on blast db download")
    parser.addoption("--travisonly", action="store_true",
        help="run travisonly tests")
    parser.addoption("--notravis", action="store_true",
        help="don't run notravis tests")


def pytest_collection_modifyitems(config, items):
    runslow =  config.getoption("--runslow")
    runconcat = config.getoption("--runconcat")
    runweb = config.getoption("--runweb")
    nolocalblast = config.getoption("--nolocalblast")
    travisonly = config.getoption("--travisonly")
    notravis = config.getoption("--notravis")

    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    skip_concat = pytest.mark.skip(reason="need --runconcat option to run")
    # skip_concat_full = pytest.mark.skip(reason="need --runconcatfull option to run")
    skip_web = pytest.mark.skip(reason="need --runweb option to run")
    skip_localblast = pytest.mark.skip(reason="--nolocalblast skips this test")
    # skip_test = pytest.mark.skip(reason="skip all for now")
    for item in items:
        if "slow" in item.keywords and not runslow:
            item.add_marker(skip_slow)
        if "concat" in item.keywords and runconcat:
            item.add_marker(skip_concat)
        # if "concatfull" in item.keywords and not runconcatfull:
        #     item.add_marker(skip_concat_full)
        if "web" in item.keywords and not runweb:
            item.add_marker(skip_web)
        #elif "fast" not in item.keywords:
        #    item.add_marker(skip_test)
        if "localblast" in item.keywords and nolocalblast:
            item.add_marker(skip_localblast)
        if "travis" in item.keywords and not travisonly:
            item.add_marker(skip_localblast)
