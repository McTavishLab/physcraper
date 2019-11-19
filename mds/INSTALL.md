[Back home](../README.md)


# Installing `physcraper`

## Preinstallation requirements

- `physcraper` currently requires python2.7.
Work is in progess on cleaning up the requirements to make it python3 ready.
- git

## Download physcraper 

```
git clone git@github.com:McTavishLab/physcraper.git
```

## Create a python virtual environment

Recommended install procedure is using a virtual environment (are there any other ways?). From your terminal do:

```
cd replace-this-by-your-own-path/physcraper
virtualenv --python=/usr/bin/python2.7 venv-physcraper
```

If that fails, you may have to use:

```
python2.7  -m pip install virtualenv
python2.7  -m virtualenv venv-physcraper
```

## Activate the installed virtual environment
Once you have a venv-physcraper directory, **_activate_** it with:

```
source venv-physcraper/bin/activate
```

Force numpy library installation. _**This is key**_:

```
pip install --force-reinstall numpy==1.14.5
pip install -r requirements.txt
python setup.py install
```

Dependencies:

Currently complete phylogenetic updating WITH `physcraper` requires
[raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) and [MUSCLE](install-muscle.md) to be installed and in the path.

You can check if they are already installed with:

```
which muscle
which raxmlHPC
```

You can generate updated alignments using other phylogenetic software by creating a Scrape object and running `scrape.align__query_seqs()`

##### Optional dependencies (need to be in path): 

- PaPaRa http://sco.h-its.org/exelixis/web/software/papara/index.html 
- BLAST+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

##### Python packages: 
These will all be installed if you install physcraper using `python setup.py install`

(but note, if you are using virtualenv there are some weird interactions with setuptools and python 2.7.6)

- Dendropy https://pythonhosted.org/DendroPy/ 
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser 

[Previous: Back home](../README.md)

[Next: Running  `physcraper`](/running.md) 