import pytest

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",
        help="run slow tests")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    # skip_test = pytest.mark.skip(reason="skip all for now")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
        #elif "fast" not in item.keywords:
        #    item.add_marker(skip_test)
