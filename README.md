# physcraper
Continual gene tree updating. 
Uses a tree from Open tree of Life and an alignment to search for and add homologous sequences to phylogenetic inference. 

![](https://cdn.rawgit.com/snacktavish/physcraper/master/docs/physcraper.svg)

Still work in progress (documentation in particular), please contact ejmctavish, gmail if you need any help!
## There is a full example python script with comments in docs/example.py 
(it takes a ehile long time though)

###Dependencies (need to be in path): 
- PaPaRa http://sco.h-its.org/exelixis/web/software/papara/index.html 
- Raxml http://sco.h-its.org/exelixis/web/software/raxml/index.html 

###Python packages: 
These will all be installed if you install physcraper using 
    python setup.py install

(but note, if you are using vritualenv there are some weird interactions with setuptools and python 2.7.6)

- Dendropy https://pythonhosted.org/DendroPy/ 
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser 

Inputs needed are:
- ott_study_id =  OpenTree study identifier  
- ott_tree_id  = Tree id from that study  
- seqaln = the sequence alignment that generated that tree  
- matrix_type = alignment matrix type (only tested with fasta so far)
- Working directory name (will be created by run)

Currently Physcraper relies on metadata information from the Open Tree of Life,
and only uses trees from that database.
Go to https://tree.opentreeoflife.org/curator to find a tree, or upload it!
You can get the tree ID by clicking on your tree of interest, and looking at the URL.

###Taxon infomation from ncbi
It is easist if you keep the taxon information in the included taxonomy folder. (the file is too big for github)
To get it from the NCBI ftp site

    rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz taxonomy/gi_taxid_nucl.dmp.gz  
    gunzip taxonomy/gi_taxid_nucl.dmp.gz




Physcraper generates an ATT (alignment, taxonomy, tree) object.
This is important because it tracks the shared namespaces across these three data objects, whcih can otherwise get a bit seperated.


