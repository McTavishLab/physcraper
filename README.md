# physcraper
Continual gene tree updating. 
Uses a tree from Open tree of Life and an alignment to search for and add homologous sequences to phylogenetic inference. 

![](https://cdn.rawgit.com/snacktavish/physcraper/master/docs/physcraper.svg)

Still work in progress (documentation in particular), please contact ejmctavish, gmail if you need any help!
## There is a full example python script with comments in docs/example.py


###Dependencies (need to be in path): 
- PaPaRa http://sco.h-its.org/exelixis/web/software/papara/index.html 
- Raxml http://sco.h-its.org/exelixis/web/software/raxml/index.html 

###Python packages: 
- Dendropy https://pythonhosted.org/DendroPy/ 
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download

Inputs needed are:
- ott_study_id =  OpenTree study identifier  
- ott_tree_id  = Tree id from that study  
- seqaln = the sequence alignment that generated that tree  
- matrix_type = alignment matrix type (only tested with fasta so far)  
- Working directory name


Currently this is relying on metadata information from Open Tree of Life,
and only uses trees from that database.
Go to https://tree.opentreeoflife.org/curator to find a tree, or upload it!
