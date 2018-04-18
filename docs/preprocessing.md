

Preprocessing steps:

Get alignement from ?treebase?, using orginal file:
convert to fasta and break out to single gene: e.g. rpb2 3218-6016 (this info is in the nexus)
using:
    preprocess.py input.nex output_stub start stop

e.g.
    
    preprocess.py M4058.nexorg rpb2 3218 6016

This will write a fasta file with just your gene/region of interest

