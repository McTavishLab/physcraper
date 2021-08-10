While Physcraper can be installed via pip,
in order to easily access the example data and ancillary files, we recommend downloading
the Physcraper repository from GitHub and installing it locally following the instructions below.
This process will also install the following python packages:

- Dendropy https://pythonhosted.org/DendroPy/
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on the Physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser


## Downloading Physcraper

First step is to download Physcraper to your computer.

You can do this with Git:

```
git clone https://github.com/McTavishLab/physcraper.git
```

or, you can download the repository from https://github.com/McTavishLab/physcraper.git

Now, move to the newly created "physcraper" directory with `cd physcraper` to continue.

Next step is to create a virtual
environment to run Physcraper on. You can do this using Anaconda
or Virtualenv.

## Anaconda virtual environment

For this option you will of course need [Anaconda](https://www.anaconda.com/products/individual) installed. You can follow installation instructions on [Anaconda's documentation
website](https://docs.anaconda.com/anaconda/install/).

Now you can create a "conda virtual environment" with:

```
   conda env create -f cond_env.yml
   conda activate physcraper_env
   pip install -r requirements.txt
   pip install -e .
```

Note the "dot" at the end of that last command, and it should be ready!

## Virtualenv virtual environment

For this option you will need [Virtualenv](https://pypi.org/project/virtualenv/) installed.

Now you can go ahead and create a "Python virtual environment".

Remember you need to be in the "physcraper" folder (go there with `cd physcraper`).
Once there do:

```
virtualenv -p python3 venv-physcraper
```

This will create a python 3 virtual environment named "venv-physcraper".

_Activate_ the virtual environment with:

```
source venv-physcraper/bin/activate
```

Finally, install Physcraper inside the virtual environment:

```
pip install -r requirements.txt
pip install -e  .
```

Do not miss the "dot" at the end of that last command!

The virtual environment remains active even if you change directories.
So, Physcraper will run from anywhere, while the virtual environment is activated.


**Note** that you will have to activate the virtual environment with `source venv-physcraper/bin/activate`
every time you want to run Physcraper.

After you are finished working with Physcraper and you don't want to run it anymore, deactivate the virtual environment with:

```
deactivate
```


## Checking for dependencies

Currently complete phylogenetic updating with Physcraper requires
[raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) and [MUSCLE](install-muscle.md) to be installed and in the path.

You can check if they are already installed with:

```
which muscle
which raxmlHPC
```

## Checking installation success on remote searches

To test a full run with pre-downloaded BLAST results, copy the example results using:

    cp -r docs/examples/pg_55_web pg_55_test

and then run:

    physcraper_run.py --study_id pg_55 --tree_id tree5864 --treebase --bootstrap_reps 10 --output pg_55_test

There is more info on all the parameter settings in the documentation section [Run](https://physcraper.readthedocs.io/en/latest/physcraper_run.html), but briefly, this gets a tree (tree5864) from study pg_55 on OpenTree, pulls the alignment from TreeBASE, blasts the sequences, and does 10 bootstrap reps on the final phylogeny.


This example tests all the components except for the actual remote BLAST searches (because they can be very slow).
To check if your installation was successful for remote searches, try running a full analysis:

    physcraper_run.py --study_id pg_55 --tree_id tree5864 --treebase --bootstrap_reps 10 --output pg_55_new

This run will take a while - once it starts blasting, that means it's working! You can use Ctrl-C to cancel.
