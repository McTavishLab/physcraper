## Details on running using physcraper_run.py


## QuickStart

For a simple run you just need the study id and tree id from OpenTree (see more about searching in FindTrees.md),
and an alignment file that goes with that tree.

    physcraper_run.py -s <study_id> -t <tree_id> -a <alignment_file_path> -as <alignment_schema> -o <output_directory>


To update this tree
https://tree.opentreeoflife.org/curator/study/view/ot_350/?tab=home&tree=Tr53296

(alignment already downloaded from treebase)


    physcraper_run.py -s ot_350 -t Tr53297 -a docs/examples/ot_350Tr53297.aln -as nexus -o ot_350_analysis


## Configuration paramaters


To see all the config paramaters, use `physcraper_run.py -h`  


### Input Data


Tree information (required)
  -s STUDY_ID, --study_id STUDY_ID
                        OpenTree study id
  -t TREE_ID, --tree_id TREE_ID
                        OpenTree tree id

Alignment information (reuired)

  -a ALIGNMENT, --alignment ALIGNMENT
                        path to alignment
  -as ALN_SCHEMA, --aln_schema ALN_SCHEMA
                        alignment schema (nexus or fasta)

OR

  -tb, --treebase       download alignment from treebase

Tree and alignment information are required.  
After an analysis has been run, they can be reloaded from a directory from a previous run.  

  -re RELOAD_FILES, --reload_files RELOAD_FILES
                        reload files and configureation from dir


REQUIRED:

  -o OUTPUT, --output OUTPUT
                        path to output directory

Optional:

  -st SEARCH_TAXON, --search_taxon SEARCH_TAXON
                        taxonomic id to constrain blast search. format ott:123
                        or ncbi:123. Deafult will use ingroup of tree on
                        OpenTree, or MRCA of input tip





### Configuration paramaters

The configuration paramaters may be set in a config file, and then passed into the analysis run. See example.config for an example.


  -c CONFIGFILE, --configfile CONFIGFILE
                        path to config file

If a config file input is combined with comand line configuration parameters, the command line values will ovverride those in the config file.


## Blast search paramaters

  -e EMAIL, --email EMAIL
                        email address for ncbi blast searches

  -r, --repeat          repeat search until no no sequences are found


  -ev EVAL, --eval EVAL
                        blast evalue cutoff
  -hl HITLIST_LEN, --hitlist_len HITLIST_LEN
                        number of blast searches to save per taxon


You can use a local blast database:
To setup see doc/LocalDB.md

  -db BLAST_DB, --blast_db BLAST_DB
                        local download of blast database




  -nt NUM_THREADS, --num_threads NUM_THREADS
                        number of threads to use in processing


You can use your own blast database, for example set up on an AWS server.

\
## Sequence filtering parameters

  -tp TRIM_PERC, --trim_perc TRIM_PERC
                        minimum percentage of seqs end of alignemnts
  -rl RELATIVE_LENGTH, --relative_length RELATIVE_LENGTH
                        max relative length of added seqs, compared to input
                        aln len
  -spn SPECIES_NUMBER, --species_number SPECIES_NUMBER
                        max number of seqs to include per species

  -de DELAY, --delay DELAY
                        how long to wait before blasting the same sequence
                        again
 
#### Tree search parameters
  -no_est, --no_estimate_tree
                        don't estimate tree (default is False)

  -bs BOOTSTRAP_REPS, --bootstrap_reps BOOTSTRAP_REPS
                        number of bootstrap reps


#### Internal arguments


  -tx TAXONOMY, --taxonomy TAXONOMY
                        path to taxonomy