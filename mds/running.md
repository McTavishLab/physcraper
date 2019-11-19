[Back home](../README.md)

# Running physcraper

You can run physcraper directly from the terminal or using a Jupyter notebook.
In any case, it will run on a python virtual environment. If you have not installed this, go to section [installing physcraper](INSTALL.md).

If you want to run `physcraper` with Jupyter notebooks, you will also need to make the virtual environment available there. This is done in different ways, depending on the python version you have installed in your virtual environment. You can verify it by doing `python --version` from your virtual environment.

## If your virtual environment has python 2

This is the only workaround I have found so far. Simply [activate your virtual environment](INSTALL.md#activate) and install jupyter notebooks within it with:

```
pip install jupyter notebook
```

## If your virtual environment has python 3

Do the following (taken from [here](https://janakiev.com/blog/jupyter-virtual-envs/) and [here](https://stackoverflow.com/questions/30604952/pip-default-behavior-conflicts-with-virtualenv)):

Activate your virtual environment

```
source venv-physcraper/bin/activate
```

Now do

```
pip install ipykernel
python -m ipykernel install --user --name=venv-physcraper
```

## If your virtual environment is in python 2

That did not work because the virtual environment is in python 2, and this is somehow set up to python 3.

The way I solved it was to 
Then, just launch it with

```
jupyter notebook
```

You can move to the directory where you are running the analysis before launching it, or you can navigate from within the notebook, once it is open.

[Previous: Installing `physcraper`](INSTALL.md)

[Next: Examples](examples.md) 
