[Back home](../README.md)

# Running physcraper

You can run `physcraper` directly from the terminal or using a Jupyter notebook.

In any case, it will run on a python virtual environment. You can find instructions on how to set this up in the section [installing physcraper](INSTALL.md).

## With Jupyter notebooks

You will have to add Jupyter notebooks to the python virtual environment. This is done in different ways, depending on the python version that has been installed in your virtual environment. You can verify this by doing `python --version` from your [activated](INSTALL.md) virtual environment.

### If your virtual environment has python 2

This is the only workaround I have found so far. Simply [activate your virtual environment](INSTALL.md) and install jupyter notebooks within it:

```
pip install jupyter notebook
```

### If your virtual environment has python 3

(taken from [here](https://janakiev.com/blog/jupyter-virtual-envs/) and [here](https://stackoverflow.com/questions/30604952/pip-default-behavior-conflicts-with-virtualenv))
Activate your virtual environment

```
source venv-physcraper/bin/activate
```

Install *ipykernel* library

```
pip install ipykernel
python -m ipykernel install --user --name=venv-physcraper
```


Now, you can launch *jupyter notebook* from the terminal

```
jupyter notebook
```

Or, open it from your applications folder and choose *venv-physcraper* as your virtual environment.




[Previous: Installing `physcraper`](INSTALL.md)

[Next: Examples](examples.md) 
