## Updating the README

The `README.Rmd` file generates the `README.md` file, which in turn generates what is shown on Physcraper's home page at
[its GitHub repository](https://github.com/McTavishLab/physcraper#readme), and
[PyPI's description](https://pypi.org/project/Physcraper/).

To update `README.md` from `README.Rmd` file, you need R and the `rmarkdown` package installed to run:

```bash
R -e 'rmarkdown::render("README.Rmd")'
```

The `index.rst` file that lives in the `docs/source/` folder controls the home page at [readthedocs](https://physcraper.readthedocs.io/en/main/index.html), which is updated automatically as you push to GitHub.

To update any of these, you have to modify, as needed, `README.Rmd` and `docs/source/index.rst`, as well as the following `.md` files living in the [docs/mds/](https://github.com/McTavishLab/physcraper/tree/main/docs/mds) folder:

- intro-badges.md
- intro-logo.md
- intro-part1.md
- citation.md
- license.md
- contact.md
- updating-the-readme.md, _aka, this file_ <span>&#10024;</span>

To create new sections, you just need to create new `.md` files in
`docs/mds/` and make sure to add them to `README.Rmd` and `docs/source/index.rst`
