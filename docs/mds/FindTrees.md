Physcraper takes as input a phylogenetic tree and the corresponding DNA alignment.
This section shows how to get one for any given taxon from the OpenTree database
using the Physcraper command line.
If you already have a tree and an alignment of your own (or downloaded from somewhere else) that you want to update
with Physcraper, please go to the [Run section - Starting with your own tree](https://physcraper.readthedocs.io/en/latest/physcraper_run.html#starting-with-your-own-tree).

## Find a tree to update from OpenTree

To search for trees on OpenTree with your taxon of interest, you can use the command
`find_trees.py`,


    find_trees.py [-h] [-t TAXON_NAME] [-ott OTT_ID] [-tb] [-o OUTPUT]

where,

  -h, --help                                Shows this list of arguments as help message and exit

  -t TAXON_NAME, --taxon_name TAXON_NAME    Specifies the name of the search taxon

  -ott OTT_ID, --ott_id OTT_ID              Specifies the OTT id number of the search taxon

  -tb, --treebase                           Returns studies with TreeBASE data only

  -o OUTPUT, --output OUTPUT                Gives the output file name and optionally the path

<br/>

For example, to find all trees in OpenTree that contain the family of flowering plants Malvaceae, do:

    find_trees.py -t Malvaceae -o malvacea.txt

If you happen to know the taxon OTT id, or you have already obtained it from the
OpenTree website [taxon homepage](https://tree.opentreeoflife.org/opentree/argus/ottol@279960/Malvaceae), you can do:

    find_trees.py -ott 279960 -o malvacea.txt


## Find the corresponding alignment on TreeBASE

To find trees with a corresponding alignment on TreeBASE use the flag `-tb` or `--treebase`:

    find_trees.py -t Malvaceae -tb -o malvacea.txt
