Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given. (This page modified from https://github.com/pyOpenSci/cookiecutter-pyopensci)

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/McTavishLab/physcraper/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

### Write Documentation

Physcraper could always use more documentation, whether as part of the
official Physcraper docs, in docstrings, or even on the web in blog posts,
articles, and such.

To edit the readthedocs documentation, modify the files in the `docs/md` folder.
You can check the documentation build locally. To do so, make sure that you are in
a `Physcraper` virtual environment and change directories to the `docs` folder with:

  cd docs

Then build the html files with:

  make html

Now you can open the html version of your documentation:

  open build/html/index.html


### Submit Feedback

The best way to send feedback is to file an issue at https://github.com/McTavishLab/physcraper/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Get Started!

Ready to contribute? Here's how to set up Physcraper for local development.

(1) Fork the Physcraper repo on GitHub https://github.com/McTavishLab/physcraper
(2) Clone your fork locally::

    $ git clone git@github.com:your_name_here/physcraper.git

(3) Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv venv-physcraper
    $ cd physcraper/
    $ pip install -e .

(4) Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

(5) When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ pytest tests

(6) Use Pylint to check your code. Move to the "bin" or "physcraper" directory to use the `.pylintrc` config file, then run::

    $ pylint insert_name_of_module_here.py

(7) Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

(8) Submit a pull request through the GitHub website.


Extra: Count the number of functions in any given module

    from inspect import getmembers, isfunction
    foos = [o for o in getmembers(physcraper) if isfunction(o[1])]
    len(foos)
