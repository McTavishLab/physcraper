[Back home](../README.md)

# Running physcraper

You can run `physcraper` directly from the terminal or using a Jupyter notebook.

In any case, it will run on a python virtual environment. If you have not installed this, go to section [installing physcraper](INSTALL.md).

If you want to run `physcraper` using Jupyter notebooks, you will have to add the python virtual environment to it. This is done in different ways, depending on the python version you have installed in your virtual environment. You can verify this by doing `python --version` from your virtual environment.

## If your virtual environment has python 2

This is the only workaround I have found so far. Simply [activate your virtual environment](INSTALL.md#L24) and install jupyter notebooks within it with:

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


Finally, you can launch jupyter notebook from the terminal with:

```
jupyter notebook
```

Or open it from your applications folder and choose venv-physcraper as your virtual environment.




[Previous: Installing `physcraper`](INSTALL.md)

[Next: Examples](examples.md) 
