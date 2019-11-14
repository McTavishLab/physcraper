# physcraper

[![Build Status](https://travis-ci.org/McTavishLab/physcraper.svg?branch=dev)](https://travis-ci.org/McTavishLab/physcraper)[![Documentation](https://readthedocs.org/projects/physcraper/badge/?version=latest&style=flat)](https://physcraper.readthedocs.io/en/latest/)[![codecov](https://codecov.io/gh/McTavishLab/physcraper/branch/dev/graph/badge.svg)](https://codecov.io/gh/McTavishLab/physcraper)

Continual gene tree updating. 
Use a tree (from the litteraturem, a synthetic tree from Open Tree of Life, or your own tree) and an alignment of any size(?) to search for and add homologous sequences for phylogenetic inference. 

![](https://cdn.rawgit.com/snacktavish/physcraper/master/docs/physcraper.svg)


The tool is under current development in the McTavish Lab.
Please contact ejmctavish, gmail if you need any help!

- [Installation](mds/INSTALL.md)
- [Run it](mds/running.md)
- [Examples]()
- [Documentation]
- [Tests]

### Requirements

see also [here](./How_to_start.md)

##### Databases

The tool uses several databases, which can automatically be downloaded/updated from ncbi. If the tool wants to access the site, you will be asked for input (`yes` or `no`), as the tool will then access ncbi, which is a US government website!


### Tutorial

For a description of which settings need to be changed and how to set-up a run, see [here](./How_to_start.md).


### Example runs and datasets

 There is a full example python script with comments in `docs/example.py`.
 Some more example files can be found in `docs/example_scripts/`.

 If you want to try running physcraper use the testdata which is in `tests/data/tiny_test_example/`

 More information will follow soon.


### Documentation

The Documentation about the different classes can be found [here](https://physcraper.readthedocs.io/en/latest/).

### Tests

There are some tests [here](./test/) and [here](./ws-test/), which test the major functionality of the code. If you want to test if the code works on your machine, please run `python tests/testfilesetup.py` and then `sh tests/run_test.sh`,  `sh ws-tests/run_ws-tests.sh`.

