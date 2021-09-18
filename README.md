
<!-- README.md is generated from README.Rmd; please edit the .Rmd file and then from R do rmarkdown::render("README.Rmd")-->

<img align="right" width="250" src="https://raw.githubusercontent.com/McTavishLab/physcraper/main/docs/physcraper-long.png">

# Welcome to Physcraper’s repository\!

[![Build
Status](https://travis-ci.org/McTavishLab/physcraper.svg?branch=main)](https://travis-ci.org/McTavishLab/physcraper)
[![Documentation](https://readthedocs.org/projects/physcraper/badge/?version=main&style=flat)](https://physcraper.readthedocs.io/en/main/)
[![codecov](https://codecov.io/gh/McTavishLab/physcraper/branch/main/graph/badge.svg)](https://codecov.io/gh/McTavishLab/physcraper)
[![pyOpenSci](https://tinyurl.com/y22nb8up)](https://github.com/pyOpenSci/software-review/issues/26)
[![DOI](https://zenodo.org/badge/41294748.svg)](https://zenodo.org/badge/latestdoi/41294748)
[![NSF-1759846](https://img.shields.io/badge/NSF-1759846-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1759846)
[![NSF-1759838](https://img.shields.io/badge/NSF-1759838-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1759838)

<p>

</p>

<p>

</p>

## Automated gene tree updating <br> with the Open Tree of Life

Use a phylogenetic tree and a DNA alignment to automatically find and
add nucleotide sequences from a genetic database, to reproducibly
improve and advance phylogenetic knowledge within a biological group.

Physcraper relies on
[taxonomic](https://tree.opentreeoflife.org/about/taxonomy-version/ott3.3)
and [phylogenetic](https://github.com/OpenTreeOfLife/phylesystem-1)
resources and [programmatic
tools](https://github.com/OpenTreeOfLife/germinator/wiki/Open-Tree-of-Life-Web-APIs)
from the [Open Tree of
Life](https://tree.opentreeoflife.org/opentree/argus/opentree12.3@ott93302)
project.

Physcraper also leverages on programmatic tools from the
[TreeBASE](https://treebase.org/treebase-web/urlAPI.html) project and
[NCBI](https://www.ncbi.nlm.nih.gov/home/develop/api/), as well as
multiple software projects listed as [requirements](#requirements)
below, to create an automatic and reproducible workflow for
phylogenetics.

You are now on the code repository. Please refer to Physcraper’s
[documentation website](https://physcraper.readthedocs.io/en/main/) for
more details
    on:

  - [installation](https://physcraper.readthedocs.io/en/main/install.html)
    instructions using an Anaconda or a Virtualenv Python virtual
    environment,
  - [function](https://physcraper.readthedocs.io/en/main/apidocs.html)
    usage,
  - finding and choosing a [starting
    dataset](https://physcraper.readthedocs.io/en/main/find-trees.html),
  - [running](https://physcraper.readthedocs.io/en/main/run.html) a full
    analysis,
  - [visualizing and
    analysing](https://physcraper.readthedocs.io/en/main/results.html)
    results,
  - real life
    [examples](https://physcraper.readthedocs.io/en/main/examples.html),
    and
  - instructions for
    [developers](https://physcraper.readthedocs.io/en/main/CONTRIBUTING.html),
    and specifically on how to [update this
    readme](https://physcraper.readthedocs.io/en/main/CONTRIBUTING.html#updating-the-readme).

:hamster: :palm\_tree: :frog: :ear\_of\_rice: :panda\_face: :tulip:
:octopus: :blossom: :whale: :mushroom: :ant: :cactus: :fish:
:maple\_leaf: :water\_buffalo: 🦠 :shell: :bug: :octocat:

## Citation

If you use Physcraper, please cite:

  - Sánchez-Reyes, L.L., M. Kandziora, & E.J McTavish. (2021).
    *Physcraper: a Python package for continually updated phylogenetic
    trees using the Open Tree of Life*. BMC Bioinformatics 22, 355. doi:
    [doi.org/10.1186/s12859-021-04274-6](https://doi.org/10.1186/s12859-021-04274-6).
    <br><br>
  - Open Tree of Life, B. Redelings, L.L. Sanchez Reyes, K.A. Cranston,
    J. Allman, M.T. Holder, & E.J. McTavish. (2019). *Open Tree of Life
    Synthetic Tree (Version 12.3)*. Zenodo. doi:
    [10.5281/zenodo.3937741](https://doi.org/10.5281/zenodo.3937741)

## License

Physcraper is made available through the [GNU General Public License
v3.0](https://github.com/McTavishLab/physcraper/blob/main/LICENSE)

## Contact

The tool is under active development in the [McTavish
Lab](https://mctavishlab.github.io/). Please post a GitHub issue
[here](https://github.com/McTavishLab/physcraper/issues) or contact
<ejmctavish@ucmerced.edu> if you need any help or have feedback.

## Updating the README

The `README.Rmd` file generates the `README.md` file, which in turn
generates what is shown on Physcraper’s home page at [its GitHub
repository](https://github.com/McTavishLab/physcraper#readme), and
[PyPI’s description](https://pypi.org/project/Physcraper/).

To update `README.md` from `README.Rmd` file, you need R and the
`rmarkdown` package installed to run:

``` bash
R -e 'rmarkdown::render("README.Rmd")'
```

The `index.rst` file that lives in the `docs/source/` folder controls
the home page at
[readthedocs](https://physcraper.readthedocs.io/en/main/index.html),
which is updated automatically as you push to GitHub.

To update any of these, you have to modify, as needed, `README.Rmd` and
`docs/source/index.rst`, as well as the following `.md` files living in
the
[docs/mds/](https://github.com/McTavishLab/physcraper/tree/main/docs/mds)
folder:

  - intro-badges.md
  - intro-logo.md
  - intro-part1.md
  - citation.md
  - license.md
  - contact.md
  - updating-the-readme.md, *aka, this file* <span>✨</span>

To create new sections, you just need to create new `.md` files in
`docs/mds/` and make sure to add them to `README.Rmd` and
`docs/source/index.rst`
