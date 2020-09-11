
## Find a tree to update
To search for trees on OpenTree with your taxon of interest, you can use
find_trees.py


usage: find_trees.py [-h] [-t TAXON_NAME] [-ott OTT_ID] [-tb] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -t TAXON_NAME, --taxon_name TAXON_NAME
                        Name of search taxon
  -ott OTT_ID, --ott_id OTT_ID
                        Name of search taxon
  -tb, --treebase       Rturn studies with treebase data only
  -o OUTPUT, --output OUTPUT
                        Output file path

e.g.

    find_trees.py --taxon_name Malvaceae --treebase -o malvacea.txt

    find_trees.py --ott_id 124219 -o orcinus.txt


## Find the corresponding alignment
