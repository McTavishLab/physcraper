[Back home](../README.md)


# I. Installing physcraper

## 1. Downloading `physcraper`

```
git clone git@github.com:McTavishLab/physcraper.git
```

## 2A. Install using conda
Install anaconda  

```
   conda env create -f cond_env.yml 
   conda activate physcraper_env
   # This next step is temprary until opentree changes are uploaded to pypi
   pip install -e git+https://github.com/OpenTreeOfLife/python-opentree@get-tree#egg=opentree

```


## 2B. Install using Virtual Env
### 1. Create a python virtual environment


```
virtualenv venv-physcraper
```


### 2. Activate the installed virtual environment

Once you have a venv-physcraper directory, **_activate_** it with:

```
source venv-physcraper/bin/activate
```
Remember that you will have to activate the virtual environment every time you want to run `physcraper`.

### 3. Install `physcraper` inside the virtual environment with

```
python setup.py install
```

This will install the following python packages also:

- Dendropy https://pythonhosted.org/DendroPy/
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser


### 4. Come out of the virtual environment:

```
deactivate
```

Do this after you are finisged working with physcraper.


# II. Checking for dependencies

Currently complete phylogenetic updating WITH `physcraper` requires
[raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) and [MUSCLE](install-muscle.md) to be installed and in the path.

You can check if they are already installed with:

```
which muscle
which raxmlHPC
```


# III. Local Databases

The tool can be run using local databases, which can be downloaded and updated from the National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov/)). 

### 1. Installing blast command line tools

To blast locally you will need to install blast command line tools first.  
Find general instructions at
https://www.ncbi.nlm.nih.gov/books/NBK279671/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/


e.g. installing blast command line tools on linux:

```
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzvf ncbi-blast-2.10.0+-x64-linux.tar.gz 
 ```
 
The binaries/scripts/executables will be installed in the `/bin` folder.

### 2. Downloading the NCBI database

If you want to download the NCBI blast database and taxonomy for faster local searches
note that the download can take several hours, depending on your internet connection.

This is what you should do:

``` 
    mkdir local_blast_db  # create the folder to save the database
    cd local_blast_db  # move to the newly created folder
    update_blastdb nt  # download the NCBI nucleotide databases
    cat *.tar.gz | tar -xvzf - -i  # unzip the nucleotide databases
    update_blastdb taxdb  # download the NCBI taxonomy database
    gunzip -cd taxdb.tar.gz | (tar xvf - )  # unzip the taxonomy database
```

### 3. Downloading the nodes and names into the physcraper/taxonomy directory

```
    cd physcraper/taxonomy
    wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' 
    gunzip -f -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)
```


[Previous: Back home](../README.md)

[Next: Running  `physcraper`](running.md)
