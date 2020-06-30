[Back home](../../README.md)

# Running physcraper

The easiest way to run physcraper is using the command line tools. This way, you can directly specify arguments, and a config file will be written down for the sake of reproducibility.

As input, you will minimally need a study and tree ids from a tree uploaded to the Open Tree of Life website:

```
physcraper_run_py -s OPENTREE_STUDY_ID -t OPENTREE_TREE_ID -o OUTPUT_DIRECTORY_NAME/AND/OR/PATH
```

If you do not specify an alignment, physcraper will try to get one of the gene alignments that generated the tree. Provide the gene alignment that you want to be updated using the `-a` command:

```
physcraper_run_py -s OPENTREE_STUDY_ID -t OPENTREE_TREE_ID -o OUTPUT/DIRECTORY/NAME/AND/OR/PATH -a PATH/TO/GENE/ALIGNMENT/NAME
```

You can also run `physcraper` from a file using python with the command:

```
python my_physcraper_run.py
```

Run it interactively with

```
python -i my_physcraper_run.py
```

There are several examples of this types of run on the next section.

Besides running `physcraper` directly from the terminal, you can also run it using a Jupyter notebook.

<!--You can find instructions on how to set a Jupyter notebook up in the section [installing physcraper](INSTALL.md).-->


Default physcraper runs use the OpenTree and NCBI web services, and are slow.
They are best run on a server or a desktop.

If you want to use a tree that is not available on OpenTree, you will need to provide a table linking the labels in your alignment and phylogeny to taxon names.

[Previous: Installation](INSTALL.md)

[Next: Example runs](examples.md)

[Go to Documentation](https://physcraper.readthedocs.io/en/latest/)
