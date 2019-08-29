
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
    fi#External data files to be kept updated
#Include some sort of path to where these live?!



rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz gi_taxid_nucl.dmp.gz
gunzip /home/ejmctavish/ncbi/gi_taxid_nucl.dmp.gz

wget http://purl.org/opentree/ott/ott2.9/ott2.9.tgz 

tar -zxvf ott2.9.tgz 


wget files.opentreeoflife.org/ott/ott3.0/ott3.0.tgz
tar -xzvf ott3.0.tgz

grep ncbi: ../ott3.0/taxonomy.tsv | sed -E -e "s/([0-9]+).+?\|.+?\|(.+?)\|.+?\|.*ncbi:([0-9]+).*/\\1,\\3,\\2/" > ott_ncbi

OR

grep ncbi: ott3.0/taxonomy.tsv | sed -r -e "s/([0-9]+).+?\|.+?\|(.+?)\|.+?\|.*ncbi:([0-9]+).*/\\1,\\3,\\2/" > physcraper/taxonomy/ott_ncbi


#need dev master version of dendropy
#also:  https://github.com/jeetsukumaran/DendroPy.git


#and requires peyotl config setup!

#and don't forget to install physcraper!
