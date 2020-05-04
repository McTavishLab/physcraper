[Back home](../README.md)


# Installing `physcraper`

## Preinstallation requirements

## Download `physcraper`

```
git clone git@github.com:McTavishLab/physcraper.git
```

# Install using conda
Install anaconda  

```
   conda env create -f cond_env.yml 
   conda activate physcraper_env
   # This next step is temprary until opentree changes are uploaded to pypi
   pip install -e git+https://github.com/OpenTreeOfLife/python-opentree@get-tree#egg=opentree

```


# INstall using Virtual Env
## Create a python virtual environment


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


# Databases

The tool can be run locally using databases, which can be downloaded and updated from the National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov/)). 

To blast locally you will need to install blast command line tools.  
Instructions at
https://www.ncbi.nlm.nih.gov/books/NBK279671/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/


e.g. on linux:
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzvf ncbi-blast-2.10.0+-x64-linux.tar.gz 
    
 The binaries are in /bin


If you want to download the blast database and taxonomy for faster local searches
NOTE: this download can take several hours, depending on your internet connection.

``` 
    mkdir local_blast_db
    update_blastdb nt
    cat *.tar.gz | tar -xvzf - -i
    update_blastdb taxdb
    gunzip -cd taxdb.tar.gz | (tar xvf - )
```

# Download the the nodes and names dowloads in tothe physcraper/taxonomy directory

```
    cd physcraper/taxonomy
    wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' 
    gunzip -f -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)
```





# Python packages:
These will all be installed if you install physcraper using `python setup.py install`


- Dendropy https://pythonhosted.org/DendroPy/
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser

## Databases


[Previous: Back home](../README.md)

[Next: Running  `physcraper`](running.md)
