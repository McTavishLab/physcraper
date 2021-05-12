Single locus analyses only provide a narrow view of the evolutionary history of a group.

After assembling individual gene data sets and phylogenies using Physcraper,
it is straigtforward to combine results (alignments and data) from those analyses to obtain species tree estimates.


The multi_loci.py script combines results from multiple single locus Physcraper runs and generates concatenated and astral input files.


Usage:

    multi_loci.py [-h] [-d MULTIPLE_RUNS_FOLDER] [-o OUTPUT] [-f {concatenate,astral}]
                     [-s {fasta,nexus}] [-m {False,True}]

Arguments:

<blockquote>
<div><dl class="option-list">
<dt><kbd><span class="option">-h </span>, <span class="option">--help </span></kbd></dt>
<dd><p>Show the help message and exit.</p>
</dd>
<dt><kbd><span class="option">-d <var>DIRECTORY_NAME</var></span>, <span class="option">--locus_runs_folder <var>DIRECTORY_NAME</var></span></kbd></dt>
<dd><p>A name (and path) of a directory containing at least two output directories from different individual Phsycraper runs.</p>
</dd>
<dt><kbd><span class="option">-o <var>DIRECTORY_NAME</var></span>, <span class="option">--output <var>DIRECTORY_NAME</var></span></kbd></dt>
<dd><p>A name (and path) for a directory to write the combined results of multiple Physcraper runs. If it exists, it will be overwritten.</p>
</dd>
<dt><kbd><span class="option">-f <var>{concatenate,astral}</var></span>, <span class="option">--format <var>{concatenate,astral}</var></span></kbd></dt>
<dd><p>Format of combined output file.</p>
</dd>
<dt><kbd><span class="option">-s <var>{fasta,nexus}</var></span>, <span class="option">--schema <var>{fasta,nexus}</var></span></kbd></dt>
<dd><p>Combined output file alignment schema.</p>
</dd>
<dt><kbd><span class="option">-m <var>{False, True}</var></span>, <span class="option">--include_missing <var>{False, True}</var></span></kbd></dt>
<dd><p>Where uneven numbers of sequences are available, concatenate with gaps.
default to <code class="docutils literal notranslate"><span class="pre">False</span></code>
.</p>
</dd>
</dl>
</div></blockquote>

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
