# Run physcraper

You can run physcraper directly from the terminal or using a Jupyter notebook.
In any case, it will run on a python virtual environment. If you have not installed it, go to section [installing physcraper]()

If you want to run physcraper with Jupyter notebooks, you will also need to make the virtual environment available there. For that, you can follow these instructions (taken from [here](https://janakiev.com/blog/jupyter-virtual-envs/) and [here](https://stackoverflow.com/questions/30604952/pip-default-behavior-conflicts-with-virtualenv)):

Activate your virtual environment

```
source venv-physcraper/bin/activate
```
Now do
```
pip install ipykernel
python -m ipykernel install --user --name=venv-physcraper
```
That actually did not work because the virtual envurinmenrt is in python 2 and this is somehow set up to python 3

The way I did it was to simply install jupyter notebooks in my virtual environment. So once vevn-physcraper has been activated, do

```
pip install jupyter notebook
```
Then, just launch it

```
jupyter notebook
```
You can move to the directory where you are running the analysis before launching it, or you can navigate from the notebook, once it is open.
