
<!-- README.md is generated from README.Rmd; please edit the .Rmd file and then from R do rmarkdown::render("README.Rmd")-->

<img align="left" width="250" src="https://raw.githubusercontent.com/McTavishLab/physcraper/main/docs/physcraper-long.png">

# Welcome to Physcraper’s repository\!

[![Build
Status](https://travis-ci.org/McTavishLab/physcraper.svg?branch=main)](https://travis-ci.org/McTavishLab/physcraper)[![Documentation](https://readthedocs.org/projects/physcraper/badge/?version=latest&style=flat)](https://physcraper.readthedocs.io/en/latest/)[![codecov](https://codecov.io/gh/McTavishLab/physcraper/branch/main/graph/badge.svg)](https://codecov.io/gh/McTavishLab/physcraper)

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
multiple software projects listed as [requirements](#Requirements)
below, to create an automatic and reproducible workflow for
phylogenetics.

You are now on the code repository. Please refer to Physcraper’s
[documentation website](https://physcraper.readthedocs.io/en/latest/)
for more details
    on:

  - [installation](https://physcraper.readthedocs.io/en/latest/install.html)
    instructions using an Anaconda or a Virtualenv Python virtual
    environment,
  - [function](https://physcraper.readthedocs.io/en/latest/apidocs.html)
    usage,
  - finding and choosing a [starting
    dataset](https://physcraper.readthedocs.io/en/latest/find-trees.html),
  - [running](https://physcraper.readthedocs.io/en/latest/run.html) a
    full analysis,
  - [visualizing and
    analysing](https://physcraper.readthedocs.io/en/latest/results.html)
    results,
  - real life
    [examples](https://physcraper.readthedocs.io/en/latest/examples.html),
    and
  - tools for
    [developers](https://physcraper.readthedocs.io/en/latest/CONTRIBUTING.html).

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

## Requirements

Physcraper requires the user to install:

  - [Anaconda](https://docs.anaconda.com/anaconda/install/) <br>
    *Anaconda Software Distribution. Computer software. Vers. 4.8.0.
    Anaconda, July. 2021. Web. <https://anaconda.com>*
  - [Virtualenv](https://pypi.org/project/virtualenv/)
  - [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) <br> *Edgar RC.
    MUSCLE: multiple sequence alignment with high accuracy and high
    throughput. Nucleic Acids Res. 2004 Mar 19;32(5):1792-7.* doi:
    [10.1093/nar/gkh340](https://doi.org/10.1093/nar/gkh340)
  - [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) <br>
    *Stamatakis, Alexandros. “RAxML version 8: a tool for phylogenetic
    analysis and post-analysis of large phylogenies.” Bioinformatics
    30.9 (2014): 1312-1313.* doi:
    [10.1093/bioinformatics/btu033](https://doi.org/10.1093/bioinformatics/btu033)
  - [BLAST
    +](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    <br> This software is only needed if using a local genetic database.
    Note that BLAST + is automatically installed when [installing
    Physcraper using
    Anaconda](https://physcraper.readthedocs.io/en/stable/install.html#anaconda-virtual-environment).
    <br> *Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+:
    architecture and applications. BMC Bioinformatics 10, 421 (2009).*
    doi:
    [10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

<br> Physcraper relies on the following Python packages that are
<b>automatically</b> installed:

  - [argparse](https://docs.python.org/3/library/argparse.html)
  - [biopython](https://biopython.org/) <br> *Cock, P. J., Antao, T.,
    Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., et al. (2009).
    Biopython: freely available Python tools for computational molecular
    biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.*
  - [configparser](https://docs.python.org/3/library/configparser.html)
  - [coverage](https://coverage.readthedocs.io/)
  - [DateTime](https://docs.python.org/3/library/datetime.html)
  - [DendroPy](https://dendropy.org/primer/index.html) <br> *Sukumaran,
    J and MT Holder. 2010. DendroPy: a Python library for phylogenetic
    computing. Bioinformatics 26: 1569-1571* doi:
    [10.1093/bioinformatics/btq228](https://doi.org/10.1093/bioinformatics/btq228)
  - [future](https://python-future.org/)
  - [m2r2](https://pypi.org/project/m2r2/)
  - [nexson](https://github.com/OpenTreeOfLife/nexson)
  - [numpy](https://numpy.org/) <br> *Harris, C. R., Millman, K. J., van
    der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., …
    Oliphant, T. E. (2020). Array programming with NumPy. Nature, 585,
    357–362.* doi:
    [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
  - [OpenTree](https://github.com/OpenTreeOfLife/python-opentree) <br>
    *Emily Jane McTavish, Luna Luisa Sanchez-Reyes, Mark T. Holder.
    (2020). OpenTree: A Python package for Accessing and Analyzing data
    from the Open Tree of Life. BioRxiv 2020.12.14.422759* doi:
    [10.1101/2020.12.14.422759](https://doi.org/10.1101/2020.12.14.422759)
  - [pandas](https://pandas.pydata.org/) <br> *McKinney, W., & others.
    (2010). Data structures for statistical computing in python. In
    Proceedings of the 9th Python in Science Conference (Vol. 445,
    pp. 51–56).*
  - [pytest](https://pytest.org/)
  - [pytest-cov](https://pytest-cov.readthedocs.io/)
  - [pytest-xdist](https://pypi.org/project/pytest-xdist/)
  - [recommonmark](https://recommonmark.readthedocs.io/)
  - [requests](https://docs.python-requests.org/) <br> *Chandra, R. V.,
    & Varanasi, B. S. (2015). Python requests essentials. Packt
    Publishing Ltd.*
  - [sh](https://amoffat.github.io/sh/)
  - [sphinx](https://www.sphinx-doc.org/)
  - [urllib3](https://urllib3.readthedocs.io/)
