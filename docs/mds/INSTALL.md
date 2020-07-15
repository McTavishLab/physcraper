## Downloading Physcraper

In order to easily access the example data and ancillary files, we recommend downloading the physcraper repo from GitHub and installing lcoally.  
First step is to clone the Physcraper repo to your computer:

```
git clone https://github.com/McTavishLab/physcraper.git
```

or dowload the repo from https://github.com/McTavishLab/physcraper.git

Then move to the newly created physcraper directory (cd physcraper) and choose a type
of installation, using conda or a python virtual environment.

### Option 1: Install Physcraper using conda

First, install [anaconda](https://www.anaconda.com/products/individual)

Then, create a conda environment

```
   conda env create -f cond_env.yml
   conda activate physcraper_env
   pip install -r requirements.txt
   pip install -e .
```

You're done with the installation with conda!

### Option 2: Install Physcraper using a python virtual environment

First, create a python virtual environment

Remeber you need to be in the physcraper folder, once there do:

```
virtualenv -p python3 venv-physcraper
```
This will create a python 3 virtual environment

Once you have a venv-physcraper directory, **_activate_** it with:

```
source venv-physcraper/bin/activate
```

You will stay in the virtual environment even if you change directories and `physcraper` should run from anywhere, while the virtual environment is activated.


**Note** that you will have to activate the virtual environment every time you want to run `physcraper`


Finally, install `physcraper` inside the virtual environment:

```
pip install -r requirements.txt
pip install -e  .
```

Note the "dot" at the end of that last command!

This will also install the following python packages:

- Dendropy https://pythonhosted.org/DendroPy/
- Peyotl https://github.com/OpenTreeOfLife/peyotl (currently needs to be on physcraper branch)
- Biopython http://biopython.org/wiki/Download
- ConfigParser


After you are finished working with physcraper and you don't want to run Physcraper anymore, deactivate the virtual environment with:

```
deactivate
```


## Checking for dependencies

Currently complete phylogenetic updating WITH `physcraper` requires
[raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) and [MUSCLE](install-muscle.md) to be installed and in the path.

You can check if they are already installed with:

```
which muscle
which raxmlHPC
```
