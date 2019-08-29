# physcraper

[![Build Status](https://travis-ci.org/McTavishLab/physcraper.svg?branch=dev)](https://travis-ci.org/McTavishLab/physcraper)[![Documentation](https://readthedocs.org/projects/physcraper/badge/?version=latest&style=flat)](https://physcraper.readthedocs.io/en/latest/)[![codecov](https://codecov.io/gh/McTavishLab/physcraper/branch/dev/graph/badge.svg)](https://codecov.io/gh/McTavishLab/physcraper)

Continual gene tree updating. 
Uses a tree from Open tree of Life (or your own tree) and an alignment to search for and adds homologous sequences to phylogenetic inference. 

![](https://cdn.rawgit.com/snacktavish/physcraper/master/docs/physcraper.svg)


The tool is under current development in the McTavish Lab.
Please contact ejmctavish, gmail if you need any help!

### Requirements

see also [here](./How_to_start.md)

##### Optional dependencies (need to be in path): 

- Muscle 
- PaPaRa http://sco.h-its.org/exelixis/web/software/papara/index.html 
- RAxML http://sco.h-its.org/exelixis/web/software/raxml/index.html 
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


### Notes to self

rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz gi_taxid_nucl.dmp.gz
gunzip /home/ejmctavish/ncbi/gi_taxid_nucl.dmp.gz

wget http://purl.org/opentree/ott/ott2.9/ott2.9.tgz 

tar -zxvf ott2.9.tgz 


wget files.opentreeoflife.org/ott/ott3.0/ott3.0.tgz
tar -xzvf ott3.0.tgz

grep ncbi: ../ott3.0/taxonomy.tsv | sed -E -e "s/([0-9]+).+?\|.+?\|(.+?)\|.+?\|.*ncbi:([0-9]+).*/\\1,\\3,\\2/" > ott_ncbi

OR

grep ncbi: ott3.0/taxonomy.tsv | sed -r -e "s/([0-9]+).+?\|.+?\|(.+?)\|.+?\|.*ncbi:([0-9]+).*/\\1,\\3,\\2/" > physcraper/taxonomy/ott_ncbi

RSYNC_COMMAND=$(s.system("rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz {}.gz".format(self.config.ncbi_dmp))
os.system("tar -xzvf {}.gz".format(self.config.ncbi_dmp)))

    if [ $? -eq 0 ]; then
        # Success do some more work!

        if [ -n "${RSYNC_COMMAND}" ]; then
            # Stuff to run, because rsync has changes
        else
            # No changes were made by rsync
        fi
    else
        # Something went wrong!
        exit 1
    fi