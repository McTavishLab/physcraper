Single locus analyses only provide a narrow view of the evolutionary history of a group.

After assembling individual gene data sets and phylogenies using Physcraper,
it is straigtforward to combine the data from those analyses to generate species tree estimates.


The multi_loci.py script can be used to combine results from single locus runs.


usage:
    multi_loci.py [-h] [-d LOCUS_RUNS_FOLDER] [-o OUTPUT] [-f {concatenate,astral}]
                     [-s {fasta,nexus}] [-m INCLUDE_MISSING]



optional arguments:
  -h, --help            show this help message and exit
  -d LOCUS_RUNS_FOLDER, --locus_runs_folder LOCUS_RUNS_FOLDER
                        folder containing results directories from individual locus runs
  -o OUTPUT, --output OUTPUT
                        folder to write to
  -f {concatenate,astral}, --format {concatenate,astral}
                        output format
  -s {fasta,nexus}, --schema {fasta,nexus}
                        putput alignment file format schema
  -m INCLUDE_MISSING, --include_missing INCLUDE_MISSING
                        Where uneven numbers of sequences are available, concatenate with gaps. 
                        default = False

           
# Astral

To generate input files for an ASTRAL species tree analysis, use -f astral

e.g.
    multi_loci.py -d tests/data/precooked/multi_loc/ -f astral -o mini_species_tree
    java -jar ../ASTRAL/astral.5.7.5.jar -i mini_species_tree/genetrees.new -a mini_species_tree/mappings.txt 






# Concatenation

To concatenate multiple loci into a single alignment use -f concatenate. 
Default settings only generate concatented loci for taxa where there is a sequence at each locus (-m False)

e.g.
    multi_loci.py -d tests/data/precooked/multi_loc/ -f concatenate -s fasta -o mini_concat


To generate concantented taxa with missing loci use -m (for missing data)

    multi_loci.py -d tests/data/precooked/multi_loc/ -f concatenate -s nexus -m -o mini_concat_gaps




# SVD quartets


    python bin/multi_loci.py -d tests/data/precooked/multi_loc/ -f svdq -m -o mini_concat2
    paup4a168_ubuntu64 mini_concat2/svdq.nex
    svdq evalq=all taxpartition=species nthreads=ncpus;
