# physcraper
Gene tree updating

python scrape_align_est.py ott_study_id ott_tree_id seqaln matrix_type runname

Arguments: 
ott_study_id =  OpenTree study identifier  
ott_tree_id  = Tree id from that study  
seqaln = the sequence alignment that generated that tree  
matrix_type = alignment matrix type (only tested with fasta so far)  
runname = a name for this run

Study needs to be in Phylesystem, get the ott_study_id and ott_tree_id from the curator app.

example:

    cd shard/ascomycota
    python /home/ejmctavish/projects/otapi/physcraper/scrape_align_est.py pg_873 tree1679 tree1679.fas fasta ascomycota


right now needs to be run in own study dir, because it generates and looks for various files.  

also, it blasts each seq, which means it is slooow and gets slower.

Dependencies
- Papara
- raxmlHPC

Python modules 
- peyotl (needs to be configured to access phylesystem)  
- Bio
- Dendropy


Preprocessing steps:

Get alignement from ?treebase?, using orginal file:
convert to fasta and break out to single gene: e.g. rpb2 3218-6016 (this info is in the nexus)
using:
    preprocess.py input.nex output_stub start stop

e.g.
    
    preprocess.py M4058.nexorg rpb2 3218 6016

This will write a fasta file with just your gene/region of interest

