[Back home](../README.md)

# Running physcraper

You can run `physcraper` directly from the terminal or using a Jupyter notebook.

You can find instructions on how to set this up in the section [installing physcraper](INSTALL.md).


The easist way to start a physcraper run is using a tree uploaded to the Open Tree of Life Project,  
and a single gene alignment for those taxa.  
By using OpenTree data tips of your tree are already ampped to taxa.

There is an example of a physcraper run based on data in OpenTree in docs/examples/data_scrape.py  

Default physcraper runs use the OpenTree and NCBI web services, and are slow. 
They are best run on a server or a desktop.


If you want to use a tree that is not available on OpenTree, you will need to provide a table linking the 
labels in your alignemnt and phylogeny to taxon names.
