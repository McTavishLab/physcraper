Single locus analyses only provide a narrow view of the evolutionary history of a group.

After assembling individual gene data sets and phylogenies using Physcraper,
it is straigtforward to combine the data from those analyses to generate species tree estimates.


The multi_loci.py script can be used to combine results from single locus runs.


usage:
    multi_loci.py [-h] [-d LOCUS_RUNS_FOLDER] [-o OUTPUT] [-f {concatenate,astral}]
                     [-s {fasta,nexus}] [-m INCLUDE_MISSING]



optional arguments:

  -h, --help            show this help message and exit
  -d MULTIPLE_RUNS_FOLDER, --locus_runs_folder MULTIPLE_RUNS_FOLDER
                        folder containing at least two output directories from individual Phsycraper runs
  -o OUTPUT, --output OUTPUT
                        folder to write to
  -f {concatenate,astral}, --format {concatenate,astral}
                        output format
  -s {fasta,nexus}, --schema {fasta,nexus}
                        putput alignment file format schema
  -m INCLUDE_MISSING, --include_missing INCLUDE_MISSING
                        Where uneven numbers of sequences are available, concatenate with gaps.
                        default = False


## Astral

To generate input files for an ASTRAL species tree analysis, (https://github.com/smirarab/ASTRAL) use -f astral.
This will generate two files in the output directory.
`genetrees.new`, a concatenation of all of the genetrees produced in individual analyses,
and `mapping.txt`, a text file linking the tip lables in each of the gene trees to taxon names.

e.g.

    multi_loci.py -d tests/data/precooked/multi_loc/ -f astral -o mini_species_tree

You can run Astral diretcly on these files
e.g.

    java -jar astral.5.7.5.jar -i mini_species_tree/genetrees.new -a mini_species_tree/mappings.txt


## Concatenation

To concatenate multiple loci into a single alignment use -f concatenate.
Default settings only generate concatenated loci for taxa where there is a sequence at each locus .

e.g.
    multi_loci.py -d tests/data/precooked/multi_loc/ -f concatenate -s fasta -o mini_concat


To generate concatenated taxa with missing loci use -m (for include missing data).

    multi_loci.py -d tests/data/precooked/multi_loc/ -f concatenate -s nexus -m -o mini_concat_gaps


This will generate a concatenated alignment in the output directory with the name 'concat.aln' in the schema selected using -s (either fasta or nexus).
Each concatenated sequences is labeled with the taxon name and an integer.

The sequences from each individual run comprising the concatenated sequence are described in "concat_info.txt" in the output directory.

## SVD quartets

To write out a concatenated Nexus file with a taxon partitions block linking sequences for the same taxa, for use in SVD quartets analyses (tutorial at http://evomics.org/learning/phylogenetics/svdquartets/) use -f svdq

This will generate a Nexus file of concatenated sequences linked together by their taxon assignment in a taxon block.
The sequences from each individual run comprising the concatenated sequence are described in "concat_info.txt" in the output directory, as above.
e.g.

    multi_loci.py -d tests/data/precooked/multi_loc/ -f svdq -m -o svdq_out

This file can be used to run SVDQ in Paup
e.g.

    paup4a168_ubuntu64 mini_concat2/svdq.nex
    svdq evalq=all taxpartition=species nthreads=ncpus;
