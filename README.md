# physcraper

[![Build Status](https://travis-ci.org/McTavishLab/physcraper.svg?branch=master)](https://travis-ci.org/McTavishLab/physcraper)[![Documentation](https://readthedocs.org/projects/physcraper/badge/?version=latest&style=flat)](https://physcraper.readthedocs.io/en/latest/)

Continual gene tree updating. 
Uses a tree from Open tree of Life (or your own tree) and an alignment to search for and adds homologous sequences to phylogenetic inference. 

![](https://cdn.rawgit.com/snacktavish/physcraper/master/docs/physcraper.svg)


The tool is under current development in the McTavish Lab.
Please contact ejmctavish, gmail if you need any help!

### Requirements

see also [here](./How_to_start.md)

##### Dependencies (need to be in path): 

- PaPaRa http://sco.h-its.org/exelixis/web/software/papara/index.html 
- Raxml http://sco.h-its.org/exelixis/web/software/raxml/index.html 
- BLAST+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

##### Python packages: 
These will all be installed if you install physcraper using 
    python setup.py install

(but note, if you are using virtualenv there are some weird interactions with setuptools and python 2.7.6)

- Dendropy https://pythonhosted.org/DendroPy/ 
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser 

##### Databases

The tool uses several databases, which can automatically be downloaded/updated from ncbi. If the tool wants to access the site, you will be asked for input (yes or no), as the tool will then access ncbi, which is a US government website!


### Tutorial

For a description of which settings need to be changed and how to set-up a run, see [here](./How_to_start.md).


### Example runs and datasets

 There is a full example python script with comments in `docs/example.py`.
 Some more example files can be found in `docs/example_scripts/`.

 If you want to try running physcraper use the testdata which is in `tests/data/tiny_test_example/`

 More information will follow soon.


### Documentation

The Documentation about the different classes can be found [here](./docs/).

### Tests

There are some tests [here](./test/) and [here](./ws-test/), which test the major functionality of the code. If you want to test if the code works on your machine, please run `python tests/testfilesetup.py` and then `sh tests/run_test.sh`,  `sh ws-tests/run_ws-tests.sh`.




