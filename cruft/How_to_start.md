
## Short introduction into what the tool actually does

`PhyScraper` is a command-line tool written in python to automatically update phylogenies. As input it needs a phylogeny, the corresponding alignment and the information about the tip names and the corresponding species names. This can either be a file provided by the user or if you have uploaded your tree to Open Tree of Life the corresponding study ID.
PhyScraper will take every input sequence and blast it against the ncbi GenBank database. Sequences that are similar to the input sequence will be added to the alignment if they are a different species and/or they are longer than existing sequences, or if they differ in at least one point mutation.
Then, `PhyScraper` will place the newly found sequences onto the tree, which is then used as a starting tree for a full RAxML run. In the next round, every newly added sequence will be blasted and this continues until no new sequence are found.
After a certain time threshold (currently 14 days), the existing sequences will be blasted again to check if new sequences can be found.
After the single-gene datasets are updated, the data can be concatenated. Either, the user specifies which sequences are combined, or the tool decides randomly which sequences to combine if there are more than a single sequence for a taxon in one of the alignments.


## Short tutorial:

### Before you can start

### 0. Install general dependencies

MAC: To install python you might require [brew](https://brew.sh/).

#### 1. Install the physcraper dependencies:

* [PaPaRa](http://sco.h-its.org/exelixis/web/software/papara/index.html) - alignment tool
* [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) - tree estimation program in a [MAC](http://www.sfu.ca/biology2/staff/dc/raxml/)
    * Make sure you do `make -f Makefile.PTHREADS.gcc` from within the RAxML folder to enable multi-core calculation
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - it's needed for filter runs and when using local BLAST databases. Setup and installation information can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK1762/).


make sure the programs are accessible from everywhere, thus add them to your PATH by going to your bash file (MAC: you can open it with `vim ~/.bash_profile`) and adding the following command line:

* UNIX: `export PATH=$PATH:/path/to/my/program`
* Windows: `set PATH=%PATH%;C:\path\to\my\program`
* MAC: `export PATH=$PATH:~/path/to/my/program` in

(! set PATH=%PATH%:  it takes the current path and sets PATH to it.)

#### 2.a) download PhyScraper using the command line:
* as a normal package: `git clone https://github.com/McTavishLab/physcraper.git`
* as a git repository: `git clone 'git@github.com:McTavishLab/physcraper.git'`

#### 2.b) install an isolated Python environment
  This is currently recommended because soem of the python packages required will be installed in older versions. You do not want to change already installed python packages on your whole computer system.
  A virtual environment will locally install the python packages as needed withput affecting anything else.

  First, install the python package `virtualenv`. You can do this form anywhere in your terminal:

  `pip install virtualenv`

  Now, generate the virtual python environment with:

  `virtualenv -p python2 name_of_your_virtual_env`  # currently it only runs on python 2.7, so you may need to just say `python` instead of `python2`, depending on your system.
  LLSR: The previous does not work for me :(

  We will call it `venv-physcraper`:

  `virtualenv --python=/usr/bin/python2.7 venv-physcraper`

  If that fails, you may have to use the following two:

      `python2.7  -m pip install virtualenv`
      `python2.7  -m virtualenv venv-physcraper`


  To actually use the virtual machine once generated, you need to activate it before doing anything else, by doing:

  `source name_of_your_virtual_env/bin/activate`

  In our case, we did:

  `source venv-physcraper/bin/activate`

  If you have finished running physcraper and you want to deactivate your environment, simply type: `deactivate`

#### 3. install python requirements and dependencies:

First, force reinstalling numpy package in a previous version:

 `pip install --force-reinstall numpy==1.14.5`

The following needs to run from within the physcraper main folder, so make sure you are in that directory:

  `cd physcraper`

Now, type:
  - `python setup.py install` This process is really fast (less than 1 min).
  - `pip install -r requirements.txt` This takes 3-5 min.

#### 4. install a local instance of the BLAST database:

<!-- decide for a BLASTing method:

Depending on the size of your tree to be updated, there are things to consider.

  * **web BLAST service**: If the tree is not too large and/or you have enough time, you can run the tool with the main settings, that uses the web BLAST service. The web service is not intended for large amounts of queries and if there are too many searchs being submitted by a user, the searches are being slowed down. Another down side is, that the species name retrieval can be difficult sometimes. Advantage is that it is the most up to date database to blast against.
  * **Amazon cloud service**: If you do not have a fast computer, there are options to pay for a pre-installed cloud service using [amazon](https://aws.amazon.com/marketplace/pp/B00N44P7L6/ref=mkt_wir_ncbi_blast).
  * **local blast database**: This is the __recommended method__, as it is the fastest and does not heavily depend on good internet connection. Especially, if the trees are bigger and/or you have a relatively fast computer, this might be the best option. Ncbi regularly publishes the databases, that can easily be downloaded and initiated. -->

  * **Install a local NCBI BLAST database:**

      General information about the BLAST database can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html.

      We are currently working on implementing a version that uses the online blast searches, but this is not yet working reliably.

      UNIX: to install the NCBI BLAST+ database follow these simple tests.

      * `open a terminal`
      * `cd /to/the/folder/of/your/future/blastdb`
      * `sudo apt-get install ncbi-blast+`
      * `wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'`  # this downloads all nt-compressed files
      * `update_blastdb nt`
      * `cat *.tar.gz | tar -xvzf - -i`  # macOS `tar` does not support the `-i` flag,  you need to use homebrew to `brew install gnu-tar` and replace the `tar` command by `gtar`
      * `blastdbcmd -db nt -info`
      * `rm *.tar.gz*`

       The last command shows you if it worked correctly. 'nt' means, we are making the nucleotide database.
       The database needs to be update regularly, the program will check the dates of your databases and will ask you to update the databases after 60 days. If your databases are older, you will be asked for input, if you want to update the databases.
       <!---
       Interactive input does not work on remote machines, to stop the program from asking, change the following line in your analysis file from `conf = ConfigObj(configfi)` to `conf = ConfigObj(configfi, interactive=False)`.
       -->
       If you want to update the databases earlier go back to step 1.

       MAC: please refer to these [up to date instructions](https://www.blaststation.com/intl/members/en/howtoblastmac.html)

       WINDOWS: Not sure yet.

  *  **Install the taxonomy database:**

      Install ncbi taxonomy database to retrieve taxon information from BLAST searches into the same directory as your blastdb from the step before.

        * `cd /to/the/folder/of/your/blastdb`
        * `wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'` # Download the taxdb archive
        * `gunzip -cd taxdb.tar.gz | (tar xvf - )`  # Install it in the BLASTDB directory
        * `rm *.tar.gz*`

  * **Install the taxonomic rank database:**
       *  `wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'`
       *  `gunzip  -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)`
       *  move files into `tests/data/`
       *  `rm taxdump.tar.gz`


  * **Update the databases:**

       The databases need to be update regularly, the program will check the dates of your databases and will ask you to update the databases after 60 days. If your databases are older, you will be asked for input, if you want to update the databases.
       Interactive input does not work on remote machines, to stop the program from asking, change the following line in your analysis file from `conf = ConfigObj(configfi)` to `conf = ConfigObj(configfi, interactive=False)`.

       If you want to update the databases earlier:

    * blast db: repeat the steps listed under 'Install a local Blast database'
    * taxonomy db: run `update_blastdb taxdb`
    * rank db: repeat the steps listed under 'install the taxonomic rank database'

### Set up a run

#### **1. edit major settings in the config file**

There is an example config file in `tests/data/localblast.config`
  * **BLAST settings**:
    * **email**:
      Please specify your email address, this is recommended/required by ncbi.
    * **e_value_thresh**:
      This is the e-value that can be retrieved from BLAST searches and is used to limit the BLAST results to sequences that are similar to the search input sequence.
      It is a parameter that describes how many hits can be expected by chance from a similar-sized database during BLAST searches. Small e-value indicate a significant match. In general, shorter sequences have lower e-values, because shorter sequences have a higher probability to occur in the database by chance. For more information please refer to the ncbi BLAST resources (e.g. \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs}).
      We used an e-value of 0.001 for all example datasets.
    * **hitlist_size**:
      This specifies the amount of sequences being returned from a BLAST search. If your phylogeny does not contain a lot of nodes, it can be a low value, which will speed-up the analysis. If the sampled lineage contains many sequences with low sequence divergence it is better to increase it to be able to retrieve all similar sequences. It is not advised to have a really low hitlist size, as this will influence the number of sequences that will be added to your alignment. Low hitlist sizes might not return all best-matches but only the first 10 even though there might be more best-matches in the database \citep{shah_misunderstood-2018}. Furthermore, for example, if the hitlist size is 10, but the phylogeny which shall be updated is sparsely sampled, this might result in an updated phylogeny, that has only the parts of the phylogeny updated, that were present in the input phylogeny. Lineages that were not present might never be added, as the 10 best hits all belong to the lineages already present.
    * **location**:
      This defines which kind of BLAST service is used. At the moment leave it as is, `local`.
       * **local**: will look for a local database under the path specified under `localblastdb`
       * **remote**: either ncbi (default) or amazon cloud service (website needs to be defined using `url_base`)
    * **localblastdb**: needed when `location = local` path to the locally stored BLAST database, needs to have a `/` at the end of the path
    * **url_base**: needed when `location = remote` and you use a Amazon cloud service
    * **num_threads**: defines how many cores shall be used for the blasting and other parallelized parts
    * **gb_id_filename**:
      If you plan to run different settings for the same phylogeny or several runs with similar phylogenies, where there might be overlap between BLAST searches, set it to true und specify the blastsubdir to be equal among runs (see last paragraph of  [section 3](#**3. More settings that can be adjusted:**)). This will share BLAST searches between runs and thus speeds up the run time of the BLAST search.
    * **delay**: define period until earlier blasted sequences are blasted again.

  * **Physcraper settings**:
    * **unmapped**: It is used when working with trees from OToL. Sometimes not all sequences of the study were mapped to OToL species identifiers.
      * **keep**: Unmapped tips will be kept and the is set to unknown_NAME_OF_THE_MRCA.
      * **remove**: Unmapped tips will be removed from the phylogeny and the alignment.
    * **seq_len_perc**: Here you can specify the minimum percentage length of newly found sequences to be added in comparison to the original alignment.
    * **trim_perc**: How many sequences need to have information at the beginning and end of an alignment, to not be trimed.
    * **max_len**: maximum length of sequences that are added to alignment, must be greater than 1.
  * **Internal settings:**
    Will follow, if you are looking for it and it's not here, send us an email.


#### **2. write your analysis file**
##### A. standard run

This is explaining how to set up a "standard run", which will add all sequences, that are similar and long enough to the alignment, as long as they are no subsequences of an already existing sequence. Optional arguments are explained in the following section.

Depending if you have uploaded your tree prior to analysis to the [OpenTreeofLife website (OToL)](https://ot14.opentreeoflife.org/opentree/argus/opentree9.1@ott93302), there are two main options:

1. using OpenTreeOfLife study ID:

      There is an example file in `docs/example_scripts/OToL_standard_run.py` it is based on the wrapper function `standard_run()`

      To obtain the study and tree ID's for an OToL run, either go to the website and query your lineage or you can run `find_studies.py` by typing in the terminal `python ./path/to/file/find_studies.py LINEAGENAME`. It will give you a studyID and a treeID, if there is a study available.

    * **study_id**: the ID of the corresponding study from OToL
    * **tree_id**: the ID of the corresponding tree from OToL
    * **seqaln**: give the path to your alignment file, must be a single gene alignment
    * **mattype**: file format of your alignment - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
    * **workdir**: path to your working directory, the folder where the intermediate and result files shall be stored.
    * **configfi**: path to your config-file, which edited in step 1.
    * **otu_jsonfi**: path to the otu json file, this will contain all the information of the sequences retrieved during the run. Usually, does not need to be edited.

2. using your own files:

    There is an example file in `docs/example_scripts/own_data_standard_local.py`, it comes with a tiny sample dataset in `tests/data/tiny_example`. The corresponding wrapper function to use in your file setup is `own_data_run()`.

    * **seqaln**: give the path to your alignment file, must be a single gene alignment
    * **mattype**: file format of your alignment - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
    * **trfn**: give the path to the file containing the corresponding phylogeny, all tips must be represented in the alignment file as well.
    * **schema_trf**: file format of your phylogeny file - currently supported: “fasta”, “newick”, “nexus”, “nexml”, “phylip”
    * **id_to_spn**: path to a comma-delimited file where tip labels correspond to species names: example file can be found in `tests/data/tiny_test_example/test_nicespl.csv`
    * **workdir**: path to your working directory, the folder where intermediate and result files shall be stored.
    * **configfi**: path to your config-file, which was edited in step 1.
    * **otu_jsonfi**: path to the otu json file, this will contain all the information of the sequences retrieved during the run. Usually, does not need to be edited.

Note: Specified paths have to start either from your root directory (e.g. `/home/USER/physcraper/path/to/file`), or can be relative from within the physcraper main folder (`./path/to/file`).

##### B. filter run:

This explains how to set up a filtered run. Here, the sequences retrieved from the BLAST search/during a standard run, are being filtered according to some user input. This is particularly of use, if there is no need to have every single sequence available being represented in your phylogeny.

Beside the standard definition, there are more input options. Currently supported are:

* **threshold**: This defines the maximum number of sequences per taxon (e.g. species) to be retrieved.

    If your input dataset already contains more sequences, there will be no additional sequences added, but also not removed.
    (If the removal of sequences that were already part of the initial phylogeny is a function someone would like to have, this should be easy to implement. Just ask.)
* **downtorank**: This defines the rank which is used to determine the maximum number of sequences per taxon.

    It can be set to None and then for all taxons, there will be the maximum number of threshold sequences retrieved.
    If it is set to species, there will no more than the maximum number of sequences randomly choosen from all sequences available for all the subspecies.
    It can be set to any ranks defined in the ncbi taxonomy browser.
* **selectby**: This defines how to select the representative sequences.

    * **blast**: All sequences belonging to a taxon will be used for a filtering blast search.

        A sequence already present in the phylogeny, or a randomly chosen sequence, will be used to blast against all other sequences from the locus with the same taxon name.
        From the sequences that pass the filtering criterium, sequences will be randomly selected as representative. The filtering criterium is that they need to be within the mean +/- standard deviation of sequence  similarity in relation to the queried sequence. See below for the explanation of the similarity value.

        If the taxon is likely monophyletic the distances will be similar and thus all sequences will fall within the mean and standard deviation of sequence similarity.
        If there are a few outlier sequences only, this seems to be likely a misidentification or mis-labeling in GenBank, outlier sequences will not be added, as they are most likely outside the allowed range of mean +/- SD.
        If the taxon is likely not monophyletic and sequences diverge a lot from each other, the mean and SD will be larger and allows to randomly pick sequences, that represent the divergence.

        As value for sequence similarity, we use bit-scores.  Bit-scores are log-scaled scores and a score is a numerical value that describes the overall quality of an alignment (thus from the blasted sequence against the other available sequences). Higher numbers correspond to higher similarity. While scores are depending on database size, the rescaled bit-scores do not. Check out https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html for more detail.
    * **length** Instead of randomly choosing between sequences that are within the criteria of the blast search using sequence divergence as a criteria, here the longest sequences will be selected.
* **blacklist**: a list of sequences, that shall not be added or were identified to be removed later. This needs to be formatted as a python list containing the GenBank identifiers (e.g. `[gi number, gi number]`). Please not, that it must be the Genbank identifiers and not the accession numbers.

1. using OpenTreeofLife study ID:
There is an example file in `docs/example_scripts/OToL_filter_run.py`. The corresponding function to use in your file setup is `filter_OTOL()`.

2. using your own files:
There is an example file in `docs/example_scripts/own_data_filter_blast.py`.  The corresponding function to use in your file setup is `filter_data_run()`.

3. Use a local folder as sequence database:

     Instead of using GenBank as the source of new sequences, we can specify a folder which contains sequences in fasta format and this folder will be used as a sequence database. Then before running a standard or filter run, sequences from that folder can be added to the alignment/phylogeny if the folder contains sequences that are similar to the sequences already present in the alignment. This is intended to be used for newly sequenced material, which is not yet published on GenBank.
     To use this you need to specify:

   * add_unpubl_seq = path to folder with the sequences
   * id_to_spn_addseq_json = path to file which translates the sequences names to species names

    There is an example file in `docs/example_scripts/own_data_localdb.py`, with an example of the input files in `tests/data/local_seqs` and `tests/data/tipnTOspn_localAdd.csv` shows how to write the comma-delimited file.


#### **3. More settings that can be adjusted:**

There are some more features that can be changed.

1. define the `ingroup_mrca`:

    Often phylogenies include outgroups, and someone might not be interested in updating that part of the tree. This can be avoided by defining the most recent common ancestor. It requires the OpenTreeOfLife identifier for the group of interest.

      You can get that ID by two different approaches:

    1. run `python scripts/get_ottid.py name_of_your_ingroup`

    2. by going to [Open Tree of Life](https://ot14.opentreeoflife.org/opentree/argus/opentree9.1@ott93302) and type in the name of the lineage and get the OTT ID at the right side of the page.

    The identifying number need to be entered here:
    1. in an OToL run: within the function  `standard_run()`/`filter_OTOL()` in your analysis file in the field for `ingroup_mrca`.

    2. in an own data run: provide ID within the function `own_data_run()`/`filter_data_run()` in your analysis file in the field for `ingroup_mrca`.

    Another aspect which needs to be considered, if your group of interest is not monophyletic and you limit the search to the mrca of the group, closely related sequences that belong for example to a different genus will not be added.
    Instead of just giving a single OToL identifier, the `ingroup_mrca` can be a list: [ottid1, ottid2, ottid3] to be able to retrieve sequences from several non-monophyletic lineages.

2. Share blast result files across runs:

    1. give the path to the folder in the wrapper function of your analysis file.

    2. in your config file: change the gb_id_filename setting to True.

    Be careful! If you have different hitlist_size defined, your blast files have different numbers of sequences saved. Sharing the folder across those different settings is not recommended!



#### **4. start to update your phylogeny:**

This should be straight forward - type in your physcraper main folder:

`python ./path/to/file/analysis-file.py`

And now you just need to wait...

#### **5. Concatenate different single-gene PhyScraper runs:**

After the single-gene PhyScraper runs were updated, the data can be combined, see for example `docs/example_scripts/concat_example.py`.
Either the program randomly decides which sequences to concatenate if there are more sequences available for one loci or the user can specify a file, which sequences shall be concatenated. An example file can be found at `tests/data/concatenation_input.csv`.

#### **6. Navigating the output:**

During a physcraper run, several files are being written out:
Here is a short introduction to what they are:

 * folder with previous_run: each Physcraper loop writes out the same sets of files, after finishing one loop, they are copied there before a new round is started

 *  all files that end with .p: are pickled files which are needed to rerun the dataset
 *  replaced_inputaln.fasta: is you input alignment, where '?' have been replaced with '-'
 *  **not_added_seq.csv**: contains newly found sequences, that passed the e-filter but where not added because of other reasons (not part of the defined mrca or to short)
 *  aln_ott.phy: used to add the newly found sequences to the alignment
 *  **physcraper.fas/physcraper.tre**: alignment and tree after updating with otuPS labels, those files can be used to relabel the tipnames, using `scripts/relabel_tree_file.py`
  *  **labelled.fas/labelled.tre**: same as physcraper.fas/tre but with different label
 * physcraper_final_notrim/trim/fas/tre: trimmed/untrimmed final dataset
 *  **taxon_sampling.csv**: list of taxon names and how often they are represented in the data
 *  **logfile**: short summary of how many sequences where added/filtered during a PhyScraper run
 *  **otu_seq_info.csv**: table with all sequences that passed evalue, and length filter and where either added or not, because they where filtered during the taxon filtering.
 * place_resolve.tre: your phylogeny with the new sequences placed onto it
 * random_resolve.tre:
 * otu_dict.json: like otu_seq_info.csv but in json format
 * **RAxML files**: files produced during a RAxML run
 * **Genbank_information_added_seq.csv**: file that contains the Genbank information of the newly added sequences

 #### **7. Common error messages:**

  * **HTTP Error 429**: Too Many Requests: you probably have too many runs running and they all try to access ncbi. Within each run there is a pause to not make ncbi angry, but if there are several runs running, the timespan between runs can get to low and then ncbi restricts you from accessing the information. Re-start and it should work.
  * **KeyError: u'taxon'**: you probably forgot to activate your virtual environment.
