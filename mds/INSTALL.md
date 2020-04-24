[Back home](../README.md)


# Installing `physcraper`

## Preinstallation requirements

## Download `physcraper`

```
git clone git@github.com:McTavishLab/physcraper.git
```

## Create a python virtual environment

Recommended install procedure is using a virtual environment (are there any other ways?). From your terminal do:

```
virtualenv venv-physcraper
```


## Activate the installed virtual environment
Once you have a venv-physcraper directory, **_activate_** it with:

```
source venv-physcraper/bin/activate
python setup.py install
```

## Dependencies

Currently complete phylogenetic updating WITH `physcraper` requires
[raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) and [MUSCLE](install-muscle.md) to be installed and in the path.

You can check if they are already installed with:

```
which muscle
which raxmlHPC
```


#### Python packages:
These will all be installed if you install physcraper using `python setup.py install`

(but note, if you are using virtualenv there are some weird interactions with setuptools and python 2.7.6)

- Dendropy https://pythonhosted.org/DendroPy/
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser

## Databases

The tool can be run locally using databases, which can be downloaded and updated from the National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov/)). 

[Previous: Back home](../README.md)

[Next: Running  `physcraper`](running.md)
