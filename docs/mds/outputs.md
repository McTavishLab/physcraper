## Output files

The analysis folder has several sub directories.
each folder is labeled with a 'tag', which by default is the alignment name, but can be set in the `physcraper_run.py` arguments.

The structure consists of:

-  inputs
    -- original tree and alignment

    -- the mapping of the labels to taxa saved as `otu_info.csv`

-  blast_run
    -- blast results for each tip in both the tree and the alignment

-  run
   -- This is where intermediate processing files, and the json formatted otu information are stored. Many fo tehse files are re-used ifthe crashes and is restarted. Make sure you use a new output directory or otherwise empty this folder if you want to change run parameters.

- outputs
   -- final tree and alignment

   -- CSV file with information about each sequence

   -- The acession numbers, taxa, and sequence lengths of matches that didn't meet the sequence length cutoffs are written out to seqlen_mismatch.txt.
