# To Update or download blast DB:

Install blast command line tools:
Instructions at
https://www.ncbi.nlm.nih.gov/books/NBK279671/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/


e.g. on linux:
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzvf ncbi-blast-2.10.0+-x64-linux.tar.gz 
    
 The binaries are in /bin


# go to your local blast dblocation
update_blastdb nt
cat *.tar.gz | tar -xvzf - -i
update_blastdb taxdb
gunzip -cd taxdb.tar.gz | (tar xvf - )




#to download taxonomy databases from ncbi, place them in the 'physcraper/taxonomy directory'
cd taxonomy
wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' 
gunzip -f -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)

