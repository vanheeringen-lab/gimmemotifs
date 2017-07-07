# GimmeMotifs

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/gimmemotifs.svg)](https://badge.fury.io/py/gimmemotifs)
[![Build Status](https://travis-ci.org/simonvh/gimmemotifs.svg?branch=master)](https://travis-ci.org/simonvh/gimmemotifs)
[![Code Health](https://landscape.io/github/simonvh/gimmemotifs/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/gimmemotifs/master)
[![Documentation Status](https://readthedocs.org/projects/gimmemotifs/badge/?version=stable)](http://gimmemotifs.readthedocs.io/en/stable/?badge=stable)

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

Now you can create a new environment for GimmeMotifs:

`$ conda create -n gimme python=3 gimmemotifs=0.11.0b0`

Before using GimmeMotifs activate the environment:

`$ source activate gimme`

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

## Frequently Asked Questions (FAQ)

#### I get the following error: "Invalid value for background argument".

Currently, this is a bug in the default configuration file. Run `gimme motifs` with the additional argument `-b gc,random`. 



