# physcraper
Gene tree updating

python scrape_align_est.py ott_study_id ott_tree_id seqaln matrix_type runname

Arguments: 
ott_study_id =  OpenTree study identifier  
ott_tree_id  = Tree id from that study  
seqaln = the sequence alignment that generated that tree  
matrix_type = alignment matrix type (only tested with fasta so far)  
runname = a name for this run


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

