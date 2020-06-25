## Details on running using physcraper_run.py


## QuickStart

For a simple run you just need the study id and tree id from OpenTree (see more about searching in FindTrees.md),
and an alignment file that goes with that tree.

    physcraper_run.py -s <study_id> -t <tree_id> -a <alignment_file_path> -as <alignment_schema> -o <output_directory>


To update this tree
https://tree.opentreeoflife.org/curator/study/view/ot_350/?tab=home&tree=Tr53296

(alignment already downloaded from treebase)


    physcraper_run.py -s ot_350 -t Tr53297 -a docs/examples/ot_350Tr53297.aln -as nexus -o ot_350_analysis