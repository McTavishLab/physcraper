##Local Databases

The BLAST tool can be run using local databases, which can be downloaded and updated from the National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov/)).

### 1. Installing BLAST command line tools

To blast locally you will need to install blast command line tools first.
Find general instructions at
https://www.ncbi.nlm.nih.gov/books/NBK279671/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/


e.g. installing BLAST command line tools on **linux**:

```
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzvf ncbi-blast-2.10.0+-x64-linux.tar.gz
 ```

The binaries/scripts/executables will be installed in the `/bin` folder.

Installing BLAST command line tools on **MAC OS** is easy, with the installer. Note, however, that the BLAST executables will be installed in `usr/local/ncbi/blast` and that you will have to add this to your path in order to be able to run the executables, by adding `export PATH=$PATH:"usr/local/ncbi/blast/bin"` to the .bash_profile

If your terminal uses zshell instead of bash, make sure you're running the .bash_profile there too.


### 2. Downloading the NCBI database

If you want to download the NCBI blast database and taxonomy for faster local searches
note that the download can take several hours, depending on your internet connection.

This is what you should do:

```
    mkdir local_blast_db  # create the folder to save the database
    cd local_blast_db  # move to the newly created folder
    update_blastdb nt  # download the NCBI nucleotide databases
    # update_blastdb.pl nt  # in MAC
    cat *.tar.gz | tar -xvzf - --ignore-zeros  # unzip the nucleotide databases
    update_blastdb taxdb  # download the NCBI taxonomy database
    # update_blastdb.pl taxdb  # in MAC
    gunzip -cd taxdb.tar.gz | (tar xvf - )  # unzip the taxonomy database
```

### 3. Downloading the nodes and names into the physcraper/taxonomy directory

```
    cd physcraper/taxonomy
    wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    gunzip -f -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)
```
