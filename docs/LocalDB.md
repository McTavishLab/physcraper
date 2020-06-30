Blast Utilities

## Create an NCBI API key



### To Update or download blast DB:
This is not necessary, but will make blast searches faster.


### Install blast command line tools:
Full instructions from ncbi at [manual](https://www.ncbi.nlm.nih.gov/books/NBK279671/) and [installation](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

On Linux :

    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzvf ncbi-blast-2.10.0+-x64-linux.tar.gz

The binaries are in /bin, and you should add them to your path


### Download

    update_blastdb nt
    cat *.tar.gz | tar -xvzf - -i
    update_blastdb taxdb
    gunzip -cd taxdb.tar.gz | (tar xvf - )




### Download taxonomy databases from ncbi, place them in the 'physcraper/taxonomy directory'

    cd taxonomy
    wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    gunzip -f -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)



## Setting up an AWS blast db

To run blast searches without NCBI's required time delays, you can set up your own server on AWS (for $).
See instructions at (AWS marketplace NCBI blast)[https://aws.amazon.com/marketplace/pp/NCBI-NCBI-BLAST/B00N44P7L6]
