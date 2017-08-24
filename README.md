# GimmeMotifs

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/gimmemotifs.svg)](https://badge.fury.io/py/gimmemotifs)
[![Build Status](https://travis-ci.org/simonvh/gimmemotifs.svg?branch=master)](https://travis-ci.org/simonvh/gimmemotifs)
[![Code Health](https://landscape.io/github/simonvh/gimmemotifs/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/gimmemotifs/master)
[![Documentation Status](https://readthedocs.org/projects/gimmemotifs/badge/?version=master)](http://gimmemotifs.readthedocs.io/en/master/?badge=master)

[![DOI](https://zenodo.org/badge/676678.svg)](https://zenodo.org/badge/latestdoi/676678)

Suite of motif tools, including a motif prediction pipeline for ChIP-seq experiments.

See [full GimmeMotifs documentation](http://gimmemotifs.readthedocs.org/) for detailed installation instructions and usage examples.

For documentation on the development version see [here](http://gimmemotifs.readthedocs.org/en/latest/).

## Easy installation

The most straightforward way to install GimmeMotifs is via [conda](https://docs.continuum.io/anaconda/) using the [bioconda](https://bioconda.github.io/) channel. If you have not used bioconda yet, first set up the necessary channels (in this order):

```
$ conda config --add channels r
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

Normally, you would be able to install GimmeMotifs with one command:

`$ conda create -n gimme gimmemotifs`

However, due to an issue with the bioconda build system, I can't release the
current stable version on bioconda. Until that is fixed, you can install it as 
follows:

```
# Create an environment called gimme with all dependencies
$ conda create -n gimme python=3 pip future numpy scipy matplotlib=2 \
statsmodels scikit-learn seaborn jinja2 bedtools pybedtools \
ucsc-genepredtobed lightning xgboost r-robustrankaggreg pillow pyyaml \
diskcache six ucsc-bigbedtobed xdg xxhash readline ghostscript homer \
gadem trawler weeder xxmotif

# Activate the environment
$ source activate gimme

# Install gimmemotifs
$ pip install https://github.com/simonvh/gimmemotifs/releases/download/0.11.1/gimmemotifs-0.11.1.tar.gz
```

Python 3 is the preferred version, however, GimmeMotifs also supports Python 2. 
Don't forget to activate the environment with `source activate gimme` whenever
you want to use GimmeMotifs.

## Quick start

### Download a genome

Create a directory to store genome files.

`$ mkdir $HOME/genomes/`

To download and index a genome (all UCSC-supported genomes):

`$ gimme genome $HOME/genomes/ hg38`

Alternatively, you can index a genome directory with chromosome FASTA files on your computer.

`$ gimme index /usr/share/genomes/hg19 hg19`

### Predict some motifs:

`$ gimme motifs my_peaks.bed -g hg38 -n my_motifs`

## Help 

* Check the [FAQ](http://gimmemotifs.readthedocs.io/en/master/faq.html#faq) for
  common issues.
* The preferred way to get support is through the Github
  [issues](https://github.com/simonvh/gimmemotifs/issues/) page
* Finally, you can reach me by [mail](simon.vanheeringen@gmail.com) or
  [Twitter](https://twitter.com/svheeringen).


