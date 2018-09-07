
## Short introduction into what the tool actually does

PhyScraper is a command-line tool written in python to automatically update phylogenies. As input it needs a phylogeny, the corresponding alignment and the information about the tip names and the corresponding species names. This can either be a file provided by the user or if you have uploaded your tree to Open Tree of Life the corresponding study ID.
PhyScraper will take every input sequence and blasts it against the ncbi GenBank database. Sequences that are similar to the input sequence will be added to the alignment, if they are a different species concept and/or they are longer than existing sequences or differ in at least one point mutation.
Then it will place the newly found sequences onto the tree, which is then used as a starting tree for a full RAxML run. In the next round, every newly added sequence will be blasted and this continues until no new sequence were found.
After a certain time threshold (currently 14 days), the existing sequences will be blasted again to check if new sequences can be found.
After the single-gene datasets are updated, the data can be concatenated. Either, the user species which sequences are combined or the tool decides randomly which sequences to combine if thre are more than a single sequence for a taxon in one of the alignments.



## Short tutorial:

### Before you can start
1. install the dependencies:
     * [PaPaRa](http://sco.h-its.org/exelixis/web/software/papara/index.html) - alignment tool
     * [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) - tree estimation program
     * [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - it's needed if you want to use the filtering approach and is recommended for the standard run as well.

2. make sure the programms are accessible from everywhere, thus add them to your PATH using the command line:
     * UNIX: `export PATH=$PATH:/path/to/my/program`
     * Windows: `set PATH=%PATH%;C:\path\to\my\program`
     * MAC: `export PATH=$PATH:~/path/to/program`

    (! set PATH=%PATH%:  it takes the current path and sets PATH to it.)

3. download PhyScraper using the command line:
     * as a normal package: `git clone https://github.com/McTavishLab/physcraper.git`
     * as a git repository: `git clone 'git@github.com:McTavishLab/physcraper.git'`

4. install python requirements and dependencies:

    run from within the physcraper main folder:

      * `python setup.py install`
      * `pip install -r requirements.txt`

5. decide for a BLASTing method:

    Depending on the size of your tree to be updated, there are things to consider.

      * **web BLAST service**: If the tree is not too large and/or you have enough time, you can run the tool with the main settings, that uses the web BLAST service. The web service is limited to 100? searches a day, for users who are above this, the searches are not terminated, but are being slowed down. Advantage is that it is the most up to date database to blast against.
      * **Amazon cloud service**: If you do not have a fast computer, there are options to pay for a pre-installed cloud service using [amazon](https://aws.amazon.com/marketplace/pp/B00N44P7L6/ref=mkt_wir_ncbi_blast).
      * **local blast database**: This is the recommended method, as it is the fastest and does not heavily depend on good internet connection. Especially, if the trees are bigger and/or you have a relatively fast computer, this might be the best option. Ncbi regularly publishes the databases, that can easily be downloaded and initiated.

        Initiating a local Blast database:

        General information about the BLAST database can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html.

        In Linux to install the BLAST database do the following (for Windows and MAC please use google to figure it out, there should be plenty of information.):

          * `open a terminal`
          * `cd /to/the/folder/of/your/future/blastdb`  
          * `sudo apt-get install ncbi-blast+` # if not already installed earlier
          * `wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'`  # this downloads all nt-compressed files
          * `update_blastdb nt`
          * `cat *.tar.gz | tar -xvzf - -i`
          * `blastdbcmd -db nt -info`
          
          The last command shows you if it worked correctly. 'nt' means, we are making the nucleotide database.
          The database needs to be update regularly, go back to step 1 as soon as there is a database update to get the most recent sequences from GenBank.

6. install ncbi taxonomy databases to circumvent crashes through internet problems during queries of the ncbi website
    
      *  install the taxonomy database
          * `cd /to/the/folder/of/your/blastdb`
          * `wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'` # Download the taxdb archive
          * `gunzip -cd taxdb.tar.gz | (tar xvf - )`  # Install it in the BLASTDB directory
    
           to update the taxonomy database later, run `perl update_blastdb.pl taxdb`
           
      * install the taxonomic rank database if you want to use the FilterBlast functions.
         *  `wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'`
         * gunzip  -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)
        
### Set up a run
1. edit major settings in the config file

    There is an example config file in `tests/data/localblast.config`

      * **email**: please specify your email address, this is recommended/required by ncbi
      * **evalue threshold**: depending on your dataset and what you want to do, it might be worse checking if the e-value threshold is how you want it.
        This is the e-value that can be retrieved from BLAST searches and is used to limit the BLAST results to sequences that are similar to the search input sequence.
      * **hitlist_size**: this is the amount of sequences being returned from a BLAST research. If your phylogeny does not contain a lot of nodes, it can be a low value. If the sampled lineage contains many sequences with low sequence divergence it is better to increase it to be able to retrieve all similar sequences.
      * **location**: this defines which kind of BLAST service is used
         * **remote**: either ncbi (default) or amazon cloud service (website needs to be defined using `url_base`)
         * **local**: will look for a local database under the path specified under localblastdb, `num_threads` defines how many cores shall be used for the blasting (the more the faster).
      * **gifilename**: if you plan to run different settings for the same phylogeny or several runs with similar phylogenies, where there might be overlap between BLAST searches, set it to true und specify the blastsubdir to be equal among runs (see next section). This will share BLAST searches between runs and thus speeds up the run time of the BLAST search.
      * **seq_len_perc**: here you can specify the minimum percentage length of newly found sequences to be added in comparison to the original alignment.

2. write your analysis file - standard run 

    This is explaining how to set up a "standard run", which will add all sequences, that are similar and long enough to the alignment, as long as they are no subsequences of an already existing sequence. Optional arguments are explained in the following section.

    Depending if you have uploaded your tree prior to analysis to the [OpenTreeofLife website (OToL)](https://ot14.opentreeoflife.org/opentree/argus/opentree9.1@ott93302), there are two main options:

    Specified paths have to start either from your root directory (e.g. `/home/USER/physcraper/path/to/file`), or can be relative from within the physcraper main folder (`./path/to/file`).

    1. using OpenTreeOfLife study ID:

        There is an example file in `docs/example.py` it is based on the wrapper function `standard_run()`

        * **study_id**: the ID of the corresponding study from OToL
        * **tree_id**: the ID of the corresponding tree from OToL
        * **seqaln**: give the path to your alignment file, must be a single gene alignment
        * **mattype**: file format of your alignment - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
        * **workdir**: path to your working directory, the folder where the intermediate and result files shall be stored.
        * **configfi**: path to your config-file, which edited in step 1.
        * **otu_jsonfi**: path to the otu json file, this will contain all the information of the sequences retrieved during the run. Usually, does not need to be edited.

    2. using your own files:

        There is an example file in `tests/tiny_standard_ownfile.py`, it comes with a tiny sample dataset in `tests/data/tiny_example`. The corresponding wrapper function to use in your file setup is `own_data_run()`.

        * **seqaln**: give the path to your alignment file, must be a single gene alignment
        * **mattype**: file format of your alignment - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
        * **trfn**: give the path to the file containing the corresponding phylogeny, all tips must be represented in the alignment file as well.
        * **schema_trf**: file format of your phylogeny file - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
        * **id_to_spn**: path to a comma-delimited file where tip labels correspond to species names: example file can be found in `tests/data/tiny_test_example/test_nicespl.csv`
        * **workdir**: path to your working directory, the folder where intermediate and result files shall be stored.
        * **configfi**: path to your config-file, which was edited in step 1.
        * **otu_jsonfi**: path to the otu json file, this will contain all the information of the sequences retrieved during the run. Usually, does not need to be edited.

3. write your analysis file - filter run:

    This explains how to set up a filtered run. Here, the sequences retrieved from the BLAST search/during a standard run, are being filtered according to some user input. This is particularly of use, if there is no need to have every single sequence available being represented in your phylogeny.

    Beside the standard definition, there are more input options. Currently supported are:

    * **threshold**: this defines the maximum number of sequences per taxonomic concept to be retrieved. If your input dataset already contains more sequences, there will be no additional sequences be added, but also not removed. (If the removal of sequences is a function someone would like to have, this should be easy to implement. Just ask.) 
    * **downtorank**: this defines the rank which is used to determine the maximum number per taxonomic concept. It can be set to None and then for all taxonomic concepts, there will be the maximum number of threshold species be retrieved. If it is set to species, there will no more than the maximum number of sequences randomly choosen from all sequences available for the subspecies. It can be set to any ranks defined in the ncbi taxonomy browser.
    * **selectby**: this defines how to select the selected sequences. For the moment only "blast" is supported, "length" is under development.
        * **blast**: All sequences belonging to a taxonomic concept will be used for a filtering blast search. Sequences will be randomly selected from a pool of sequences, that met the criteria of being within the mean +/- standard deviation of sequence  similarity in relation to the queried sequence (a sequence of the same taxonomic concept which is already part of the phylogeny). If the taxonomic concept is likely monophyletic the distances will be similar and thus all sequences will fall within the mean and standard deviation of sequence similarity. If there are a few outlier sequences only, this seems to be likely a misidentification or mis-labeling in GenBank, outlier sequences will not be added, as they are most likely outside the allowed range of mean +/- SD. If the taxon is likely not monophyletic and sequences diverge a lot from each other, the mean and sd value will be larger and allows to randomly pick sequences, that represent the divergence.
        * **length** Instead of randomly choosing between sequences that are within the criteria of the blast search using sequence divergence as a criteria, here the longest sequences will be selected.
    * **blacklist**: a list of sequences, that shall not be added or were identified to be removed later. This needs to be formatted as a python list containing the GenBank identifiers (e.g. `[gi number, gi number]`). Please not, that it must be the Genbank identifiers and not the accession numbers.

    1. using OpenTreeofLife study ID:
    There is an example file in `tests/filter_OTOL.py`. The corresponding function to use in your file setup is `filter_OTOL()`.

     2. using your own files:
     There is an example file in `tests/tiny_filter_ownfile.py`.  The corresponding function to use in your file setup is `filter_data_run()`.

4. start to update your phylogeny:

    This should be straight forward - type in your physcraper main folder:

    `python ./path/to/file/analysis-file.py`

    And now you just need to wait...

5. more hidden features that can be changed:

    There are some more features that can be changed if you know where, we will change the code hopefully soon, to make that easier adjustable.

    * time lapse for blasting: at the moment this is set to be 14 days. If you want to adjust the timing change `run_blast()` in the wrapper to `run_blast(delay = your_value)`
    * trim method: by default sequences will be trimmed from the alignment if it has not at least 75% of the total sequence length. This can be changed in `./physcraper/__init__.py`, in the function `trim()` the value for `taxon_missingness`. 
    * change the most recent common ancestor (mrca): often phylogenies include outgroups, and someone might not be interested in updating that part of the tree. This can be avoided by defining the most recent common ancestor. It requires the OpenTreeOfLife identifier for the group of interest. 
        
        You can get that ID by two different approaches:
        1. run `python scripts/get_ottid.py name_of_your_ingroup`
        2. by going to [Open Tree of Life](https://ot14.opentreeoflife.org/opentree/argus/opentree9.1@ott93302) and type in the name of the lineage and get the OTT ID at the right side of the page. That number needs to be provided in the corresponding wrapper function, as following:
        
        The identifying number need to be entered here:
        1. in an OToL run: within the function  `standard_run()`/`filter_OTOL()` in your analysis file in the field for `ingroup_mrca`.

        2. in an own data run: provide ID within the function `own_data_run()`/`filter_data_run()` in your analysis file in the field for `ingroup_mrca`.

        Another aspect which needs to be considered, if your group of interest is not monophyletic and you limit the search to the mrca of the group, closely related sequences that belong for example to a different genus will not be added.
    * sharing blast result files across runs: 
 
        1. give the path to the folder in the wrapper function of your analysis file.
        2. in your config file: change the gifilename setting to True. 
        
        Be careful! If you have different hitlist_size defined, your blast files have different numbers of sequences saved. Sharing the folder across those different settings is not recommended!

6. Concatenate different single-gene PhyScraper runs:
    
    After the single-gene PhyScraper runs were updated, the data can be combined, see for example `tests/data/concat_runs.py`.
    Either the program randomly decides which sequences to concatenate if there are more sequences available for one loci. Or the user can specify a file, which sequences shall be concatenated. An example file can be found at `tests/data/concatenation_input.csv`.